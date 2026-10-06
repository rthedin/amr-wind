#include "src/equation_systems/icns/source_terms/MetMastForcing.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_REAL.H"
#include "AMReX_Utility.H"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

namespace {

template <typename T>
void copy_to_device(
    const amrex::Vector<T>& host, amrex::Gpu::DeviceVector<T>& device)
{
    device.resize(host.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
}

} // namespace

// Named namespace: nvcc rejects device lambdas in functions with internal
// linkage
namespace metmast {

//! Device view of the station profiles
struct StationView
{
    const amrex::Real* x;
    const amrex::Real* y;
    const int* offset;
    const amrex::Real* z;
    const amrex::Real* vel;
    const amrex::Real* sigma;
    int nstations;
};

struct ForcingParams
{
    amrex::Real inv_rh2;
    //! Horizontal weight of the body-force footprint average
    amrex::Real inv_rh2_avg;
    amrex::Real inv_rz2;
    amrex::Real cutoff;
    amrex::Real sigma_factor;
    amrex::Real inv_tau;
    amrex::Real max_rate;
    //! Half the range-gate length of the footprint average (0: not used)
    amrex::Real half_gate;
};

/** Add the met-mast forcing to the momentum source
 *
 *  With terrain, the height is measured from the local terrain and the cells
 *  inside the terrain are not forced.
 */
template <bool HasTerrain>
void add_forcing(
    const amrex::Geometry& geom,
    const StationView& st,
    const ForcingParams& prm,
    const amrex::MultiArray4<amrex::Real const>& vel_arrs,
    const amrex::MultiArray4<amrex::Real const>& terrain_arrs,
    const amrex::MultiArray4<int const>& blank_arrs,
    amrex::MultiFab& src_term)
{
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            if constexpr (HasTerrain) {
                if (blank_arrs[nbx](i, j, k) == 1) {
                    return;
                }
                z = amrex::max<amrex::Real>(
                    z - terrain_arrs[nbx](i, j, k), 0.5_rt * dx[2]);
            }
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const auto& vel = vel_arrs[nbx];

            amrex::Real sum_w = 0.0_rt;
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> sum_wt{
                AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)};
            for (int s = 0; s < st.nstations; ++s) {
                const amrex::Real xs = x - st.x[s];
                const amrex::Real ys = y - st.y[s];
                const amrex::Real rh2 = ((xs * xs) + (ys * ys)) * prm.inv_rh2;
                if (rh2 > prm.cutoff) {
                    continue;
                }
                // Bracket the height in the station profile; outside the
                // measured range use the end value and taper the weight
                const int lbot = st.offset[s];
                const int ltop = st.offset[s + 1] - 1;
                int la = lbot;
                int lb = lbot;
                amrex::Real frac = 0.0_rt;
                amrex::Real dz_out = 0.0_rt;
                if (z <= st.z[lbot]) {
                    dz_out = st.z[lbot] - z;
                } else if (z >= st.z[ltop]) {
                    la = ltop;
                    lb = ltop;
                    dz_out = z - st.z[ltop];
                } else {
                    lb = lbot + 1;
                    while (st.z[lb] < z) {
                        ++lb;
                    }
                    la = lb - 1;
                    frac = (z - st.z[la]) / (st.z[lb] - st.z[la]);
                }
                const amrex::Real ri2 = rh2 + (dz_out * dz_out * prm.inv_rz2);
                if (ri2 > prm.cutoff) {
                    continue;
                }
                const amrex::Real w = std::exp(-0.25_rt * ri2);
                sum_w += w;
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    const int ia = (la * AMREX_SPACEDIM) + n;
                    const int ib = (lb * AMREX_SPACEDIM) + n;
                    const amrex::Real ref =
                        ((1.0_rt - frac) * st.vel[ia]) + (frac * st.vel[ib]);
                    const amrex::Real band =
                        prm.sigma_factor * (((1.0_rt - frac) * st.sigma[ia]) +
                                            (frac * st.sigma[ib]));
                    const amrex::Real target = amrex::min(
                        amrex::max(vel(i, j, k, n), ref - band), ref + band);
                    sum_wt[n] += w * target;
                }
            }
            if (sum_w <= 0.0_rt) {
                return;
            }
            const amrex::Real rate = amrex::min(
                amrex::min(sum_w, 1.0_rt) * prm.inv_tau, prm.max_rate);
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                src_arrs[nbx](i, j, k, n) -=
                    rate * (vel(i, j, k, n) - (sum_wt[n] / sum_w));
            }
        });
}

/** Weight of station s at (x, y, z) and the station levels bracketing z
 *
 *  Same weights as add_forcing, with the horizontal scale inv_rh2: the levels
 *  la and lb carry (1 - frac) and frac of the weight, and outside the measured
 *  range the end level carries all of it with a Gaussian taper.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real station_weight(
    const StationView& st,
    const ForcingParams& prm,
    const amrex::Real inv_rh2,
    const int s,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real z,
    int& la,
    int& lb,
    amrex::Real& frac)
{
    const amrex::Real xs = x - st.x[s];
    const amrex::Real ys = y - st.y[s];
    const amrex::Real rh2 = ((xs * xs) + (ys * ys)) * inv_rh2;
    if (rh2 > prm.cutoff) {
        return 0.0_rt;
    }
    const int lbot = st.offset[s];
    const int ltop = st.offset[s + 1] - 1;
    la = lbot;
    lb = lbot;
    frac = 0.0_rt;
    amrex::Real dz_out = 0.0_rt;
    if (z <= st.z[lbot]) {
        dz_out = st.z[lbot] - z;
    } else if (z >= st.z[ltop]) {
        la = ltop;
        lb = ltop;
        dz_out = z - st.z[ltop];
    } else {
        lb = lbot + 1;
        while (st.z[lb] < z) {
            ++lb;
        }
        la = lb - 1;
        frac = (z - st.z[la]) / (st.z[lb] - st.z[la]);
    }
    const amrex::Real ri2 = rh2 + (dz_out * dz_out * prm.inv_rz2);
    if (ri2 > prm.cutoff) {
        return 0.0_rt;
    }
    return std::exp(-0.25_rt * ri2);
}

/** Accumulate the footprint sums of every station level on one AMR level
 *
 *  For each level: sum of weights, weighted velocity and weighted squared
 *  velocity, 1 + 2 AMREX_SPACEDIM values. Cells covered by a finer level and
 *  cells inside the terrain are skipped.
 */
template <bool HasTerrain>
void accumulate_footprint(
    const amrex::Geometry& geom,
    const StationView& st,
    const ForcingParams& prm,
    const amrex::MultiFab& vel_mf,
    const amrex::iMultiFab& level_mask,
    const amrex::MultiArray4<amrex::Real const>& terrain_arrs,
    const amrex::MultiArray4<int const>& blank_arrs,
    amrex::Real* sums)
{
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const amrex::Real vol = dx[0] * dx[1] * dx[2];
    auto const& vel_arrs = vel_mf.const_arrays();
    auto const& mask_arrs = level_mask.const_arrays();

    amrex::ParallelFor(
        vel_mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            if (mask_arrs[nbx](i, j, k) == 0) {
                return;
            }
            amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            if constexpr (HasTerrain) {
                if (blank_arrs[nbx](i, j, k) == 1) {
                    return;
                }
                z = amrex::max<amrex::Real>(
                    z - terrain_arrs[nbx](i, j, k), 0.5_rt * dx[2]);
            }
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const auto& vel = vel_arrs[nbx];
            const int nsum = 1 + (2 * AMREX_SPACEDIM);

            for (int s = 0; s < st.nstations; ++s) {
                int la = 0;
                int lb = 0;
                amrex::Real frac = 0.0_rt;
                const amrex::Real w = station_weight(
                    st, prm, prm.inv_rh2_avg, s, x, y, z, la, lb, frac);
                if (w <= 0.0_rt) {
                    continue;
                }
                const amrex::GpuArray<int, 2> lvl{la, lb};
                const amrex::GpuArray<amrex::Real, 2> wl{
                    (1.0_rt - frac) * w * vol, frac * w * vol};
                for (int m = 0; m < 2; ++m) {
                    if (wl[m] <= 0.0_rt) {
                        continue;
                    }
                    const int base = lvl[m] * nsum;
                    amrex::Gpu::Atomic::AddNoRet(&sums[base], wl[m]);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        const amrex::Real u = vel(i, j, k, n);
                        amrex::Gpu::Atomic::AddNoRet(
                            &sums[base + 1 + n], wl[m] * u);
                        amrex::Gpu::Atomic::AddNoRet(
                            &sums[base + 1 + AMREX_SPACEDIM + n],
                            wl[m] * u * u);
                    }
                }
            }
        });
}

/** Accumulate range-gate footprint sums, like a lidar measures
 *
 *  Each gate averages the cells within half_gate of its height, with the
 *  horizontal averaging weight, and nothing beyond the end gates. The sums
 *  have the same layout as accumulate_footprint.
 */
template <bool HasTerrain>
void accumulate_gates(
    const amrex::Geometry& geom,
    const StationView& st,
    const ForcingParams& prm,
    const amrex::MultiFab& vel_mf,
    const amrex::iMultiFab& level_mask,
    const amrex::MultiArray4<amrex::Real const>& terrain_arrs,
    const amrex::MultiArray4<int const>& blank_arrs,
    amrex::Real* sums)
{
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const amrex::Real vol = dx[0] * dx[1] * dx[2];
    auto const& vel_arrs = vel_mf.const_arrays();
    auto const& mask_arrs = level_mask.const_arrays();

    amrex::ParallelFor(
        vel_mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            if (mask_arrs[nbx](i, j, k) == 0) {
                return;
            }
            amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            if constexpr (HasTerrain) {
                if (blank_arrs[nbx](i, j, k) == 1) {
                    return;
                }
                z = amrex::max<amrex::Real>(
                    z - terrain_arrs[nbx](i, j, k), 0.5_rt * dx[2]);
            }
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const auto& vel = vel_arrs[nbx];
            const int nsum = 1 + (2 * AMREX_SPACEDIM);

            for (int s = 0; s < st.nstations; ++s) {
                const amrex::Real xs = x - st.x[s];
                const amrex::Real ys = y - st.y[s];
                const amrex::Real rh2 =
                    ((xs * xs) + (ys * ys)) * prm.inv_rh2_avg;
                if (rh2 > prm.cutoff) {
                    continue;
                }
                const amrex::Real wv = std::exp(-0.25_rt * rh2) * vol;
                for (int l = st.offset[s]; l < st.offset[s + 1]; ++l) {
                    if (std::abs(z - st.z[l]) > prm.half_gate) {
                        continue;
                    }
                    const int base = l * nsum;
                    amrex::Gpu::Atomic::AddNoRet(&sums[base], wv);
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        const amrex::Real u = vel(i, j, k, n);
                        amrex::Gpu::Atomic::AddNoRet(
                            &sums[base + 1 + n], wv * u);
                        amrex::Gpu::Atomic::AddNoRet(
                            &sums[base + 1 + AMREX_SPACEDIM + n], wv * u * u);
                    }
                }
            }
        });
}

/** Add the body-force body force to the momentum source
 *
 *  The force of each station level is interpolated to the cell with the
 *  station weights; it does not depend on the local velocity.
 */
template <bool HasTerrain>
void add_body_force(
    const amrex::Geometry& geom,
    const StationView& st,
    const ForcingParams& prm,
    const amrex::Real* force,
    const amrex::MultiArray4<amrex::Real const>& terrain_arrs,
    const amrex::MultiArray4<int const>& blank_arrs,
    amrex::MultiFab& src_term)
{
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            if constexpr (HasTerrain) {
                if (blank_arrs[nbx](i, j, k) == 1) {
                    return;
                }
                z = amrex::max<amrex::Real>(
                    z - terrain_arrs[nbx](i, j, k), 0.5_rt * dx[2]);
            }
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);

            amrex::Real sum_w = 0.0_rt;
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> sum_wf{
                AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)};
            for (int s = 0; s < st.nstations; ++s) {
                int la = 0;
                int lb = 0;
                amrex::Real frac = 0.0_rt;
                const amrex::Real w = station_weight(
                    st, prm, prm.inv_rh2, s, x, y, z, la, lb, frac);
                if (w <= 0.0_rt) {
                    continue;
                }
                sum_w += w;
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    sum_wf[n] +=
                        w *
                        (((1.0_rt - frac) * force[(la * AMREX_SPACEDIM) + n]) +
                         (frac * force[(lb * AMREX_SPACEDIM) + n]));
                }
            }
            if (sum_w <= 0.0_rt) {
                return;
            }
            const amrex::Real norm = 1.0_rt / amrex::max(sum_w, 1.0_rt);
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                src_arrs[nbx](i, j, k, n) += sum_wf[n] * norm;
            }
        });
}

} // namespace metmast

MetMastForcing::MetMastForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_sim(sim)
{
    amrex::ParmParse pp_abl("ABL");
    std::string point_file;
    amrex::Vector<std::string> profile_files;
    pp_abl.query("metmast_1dprofile_file", point_file);
    pp_abl.queryarr("metmast_profile_files", profile_files);
    if (point_file.empty() && profile_files.empty()) {
        amrex::Abort(
            "MetMastForcing: set ABL.metmast_1dprofile_file and/or "
            "ABL.metmast_profile_files");
    }
    if (!point_file.empty()) {
        read_point_file(point_file);
    }
    for (const auto& fname : profile_files) {
        read_profile_file(fname);
    }
    if (num_stations() == 0) {
        amrex::Abort(
            "MetMastForcing: no measurements found in the input files");
    }

    pp_abl.query("meso_timescale", m_timescale);
    pp_abl.query("metmast_timescale", m_timescale);
    pp_abl.query("metmast_horizontal_radius", m_horizontal_radius);
    pp_abl.query("metmast_vertical_radius", m_vertical_radius);
    pp_abl.query("metmast_damping_radius", m_damping_radius);
    pp_abl.query("metmast_sigma_factor", m_sigma_factor);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_timescale > 0.0_rt, "MetMastForcing: time scale must be positive");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (m_horizontal_radius > 0.0_rt) && (m_vertical_radius > 0.0_rt),
        "MetMastForcing: radii must be positive");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_sigma_factor >= 0.0_rt,
        "MetMastForcing: sigma factor must not be negative");

    std::string forcing_type{"relaxation"};
    pp_abl.query("metmast_forcing_type", forcing_type);
    if (forcing_type == "body_force") {
        m_body_force_mode = true;
    } else if (forcing_type == "monitor") {
        m_body_force_mode = true;
        m_monitor_only = true;
    } else if (forcing_type != "relaxation") {
        amrex::Abort(
            "MetMastForcing: ABL.metmast_forcing_type must be relaxation, "
            "body_force or monitor");
    }
    pp_abl.query("metmast_averaging_time", m_averaging_time);
    m_averaging_radius = m_horizontal_radius;
    pp_abl.query("metmast_averaging_radius", m_averaging_radius);
    pp_abl.query("metmast_gate_length", m_gate_length);
    pp_abl.query("metmast_gain_schedule", m_gain_schedule);
    pp_abl.query("metmast_integral_factor", m_integral_factor);
    pp_abl.query("metmast_integral_ratio", m_integral_ratio);
    pp_abl.query("metmast_min_speed", m_min_speed);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (m_integral_factor > 0.0_rt) && (m_integral_ratio > 0.0_rt) &&
            (m_min_speed > 0.0_rt),
        "MetMastForcing: gain schedule factors and minimum speed must be "
        "positive");
    m_integral_timescale = 4.0_rt * m_timescale;
    pp_abl.query("metmast_integral_timescale", m_integral_timescale);
    pp_abl.query("metmast_force_vertical", m_force_vertical);
    pp_abl.query("metmast_max_force", m_max_force);
    pp_abl.query("metmast_output_frequency", m_output_frequency);
    pp_abl.query("metmast_restart_state", m_restart_state);
    {
        amrex::ParmParse pp_io("io");
        pp_io.query("post_processing_directory", m_post_dir);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (m_averaging_time > 0.0_rt) && (m_averaging_radius > 0.0_rt),
        "MetMastForcing: averaging time and radius must be positive");

    const int nvals = num_levels() * AMREX_SPACEDIM;
    m_footprint_mean.resize(nvals, 0.0_rt);
    m_second_moment.resize(nvals, 0.0_rt);
    m_integral.resize(nvals, 0.0_rt);
    m_force.resize(nvals, 0.0_rt);
    m_level_tau.resize(num_levels(), m_timescale);
    m_level_ti.resize(num_levels(), m_integral_timescale);
    copy_to_device(m_force, m_force_d);
    m_sums_d.resize(
        static_cast<amrex::Long>(num_levels()) * (1 + (2 * AMREX_SPACEDIM)));

    copy_to_device(m_station_x, m_station_x_d);
    copy_to_device(m_station_y, m_station_y_d);
    copy_to_device(m_level_offset, m_level_offset_d);
    copy_to_device(m_level_z, m_level_z_d);
    copy_to_device(m_level_vel, m_level_vel_d);
    copy_to_device(m_level_sigma, m_level_sigma_d);
}

MetMastForcing::~MetMastForcing() = default;

void MetMastForcing::add_level(
    const amrex::Real z,
    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& vel,
    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& sigma)
{
    m_level_z.push_back(z);
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        m_level_vel.push_back(vel[n]);
        m_level_sigma.push_back(sigma[n]);
    }
}

void MetMastForcing::read_point_file(const std::string& fname)
{
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("MetMastForcing: cannot find met-mast file " + fname);
    }
    //! x y z u v w T; the temperature is not used
    amrex::Real x, y, z, u, v, w, temp;
    while (ifh >> x >> y >> z >> u >> v >> w >> temp) {
        m_station_x.push_back(x);
        m_station_y.push_back(y);
        add_level(z, {u, v, w}, {0.0_rt, 0.0_rt, 0.0_rt});
        m_level_offset.push_back(static_cast<int>(m_level_z.size()));
    }
    if (!ifh.eof()) {
        amrex::Abort("MetMastForcing: cannot parse met-mast file " + fname);
    }
}

void MetMastForcing::read_profile_file(const std::string& fname)
{
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("MetMastForcing: cannot find profile file " + fname);
    }
    amrex::Real x, y;
    if (!(ifh >> x >> y)) {
        amrex::Abort("MetMastForcing: cannot read location in " + fname);
    }
    //! z u v w su sv sw, with z above the terrain in ascending order
    const auto nlevels_old = m_level_z.size();
    amrex::Real z, u, v, w, su, sv, sw;
    while (ifh >> z >> u >> v >> w >> su >> sv >> sw) {
        if ((m_level_z.size() > nlevels_old) && (z <= m_level_z.back())) {
            amrex::Abort(
                "MetMastForcing: heights must be increasing in " + fname);
        }
        if ((su < 0.0_rt) || (sv < 0.0_rt) || (sw < 0.0_rt)) {
            amrex::Abort(
                "MetMastForcing: negative standard deviation in " + fname);
        }
        add_level(z, {u, v, w}, {su, sv, sw});
    }
    if (!ifh.eof()) {
        amrex::Abort("MetMastForcing: cannot parse profile file " + fname);
    }
    if (m_level_z.size() == nlevels_old) {
        amrex::Abort("MetMastForcing: no heights in profile file " + fname);
    }
    m_station_x.push_back(x);
    m_station_y.push_back(y);
    m_level_offset.push_back(static_cast<int>(m_level_z.size()));
}

amrex::Real MetMastForcing::footprint_sigma(const int ilev, const int n) const
{
    const int idx = (ilev * AMREX_SPACEDIM) + n;
    return std::sqrt(
        amrex::max(
            m_second_moment[idx] -
                (m_footprint_mean[idx] * m_footprint_mean[idx]),
            0.0_rt));
}

void MetMastForcing::update_body_force(
    const metmast::StationView& st, const metmast::ForcingParams& prm) const
{
    BL_PROFILE("kynema-sgf::MetMastForcing::update_body_force");
    if (!m_initialized && !m_restart_state.empty()) {
        read_state(m_restart_state);
    }

    // Footprint sums over all AMR levels, finest level only where refined
    const int nsum = 1 + (2 * AMREX_SPACEDIM);
    const int nlevels = num_levels();
    amrex::Vector<amrex::Real> sums(
        static_cast<amrex::Long>(nlevels) * nsum, 0.0_rt);
    copy_to_device(sums, m_sums_d);

    const auto& repo = m_sim.repo();
    const bool has_terrain = repo.field_exists("terrain_height");
    const auto& velocity = m_velocity.state(FieldState::Old);
    const int finest_level = m_mesh.finestLevel();
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = amrex::makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), m_mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }
        const auto& geom = m_mesh.Geom(lev);
        const auto terrain_arrs =
            has_terrain ? repo.get_field("terrain_height")(lev).const_arrays()
                        : amrex::MultiArray4<amrex::Real const>();
        const auto blank_arrs =
            has_terrain
                ? repo.get_int_field("terrain_blank")(lev).const_arrays()
                : amrex::MultiArray4<int const>();
        if (m_gate_length > 0.0_rt) {
            if (has_terrain) {
                metmast::accumulate_gates<true>(
                    geom, st, prm, velocity(lev), level_mask, terrain_arrs,
                    blank_arrs, m_sums_d.data());
            } else {
                metmast::accumulate_gates<false>(
                    geom, st, prm, velocity(lev), level_mask, terrain_arrs,
                    blank_arrs, m_sums_d.data());
            }
        } else if (has_terrain) {
            metmast::accumulate_footprint<true>(
                geom, st, prm, velocity(lev), level_mask, terrain_arrs,
                blank_arrs, m_sums_d.data());
        } else {
            metmast::accumulate_footprint<false>(
                geom, st, prm, velocity(lev), level_mask, terrain_arrs,
                blank_arrs, m_sums_d.data());
        }
    }
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, m_sums_d.begin(), m_sums_d.end(),
        sums.begin());
    amrex::ParallelDescriptor::ReduceRealSum(
        sums.data(), static_cast<int>(sums.size()));

    // Time filter and proportional-integral controller per level
    const amrex::Real dt = m_time.delta_t();
    const bool advance = m_initialized && (dt > 0.0_rt);
    const amrex::Real alpha =
        advance ? amrex::min(dt / m_averaging_time, 1.0_rt) : 1.0_rt;
    const amrex::Real inv_ti = (m_integral_timescale > 0.0_rt)
                                   ? 1.0_rt / m_integral_timescale
                                   : 0.0_rt;
    for (int l = 0; l < nlevels; ++l) {
        const amrex::Real wsum = sums[static_cast<amrex::Long>(l) * nsum];
        if (wsum <= 0.0_rt) {
            // No cells in this footprint, e.g. inside the terrain
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                m_force[(l * AMREX_SPACEDIM) + n] = 0.0_rt;
            }
            continue;
        }
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            const int idx = (l * AMREX_SPACEDIM) + n;
            const amrex::Real m1 = sums[(l * nsum) + 1 + n] / wsum;
            const amrex::Real m2 =
                sums[(l * nsum) + 1 + AMREX_SPACEDIM + n] / wsum;
            m_footprint_mean[idx] += alpha * (m1 - m_footprint_mean[idx]);
            m_second_moment[idx] += alpha * (m2 - m_second_moment[idx]);
        }

        // Gains of this level: fixed, or scheduled from the time the air
        // spends under the force, T_p = R_h / U, which is also the delay
        amrex::Real inv_tau_l = prm.inv_tau;
        amrex::Real inv_ti_l = inv_ti;
        if (m_gain_schedule) {
            const int i0 = l * AMREX_SPACEDIM;
            const amrex::Real ux = m_footprint_mean[i0];
            const amrex::Real uy = m_footprint_mean[i0 + 1];
            const amrex::Real speed =
                amrex::max(std::sqrt((ux * ux) + (uy * uy)), m_min_speed);
            const amrex::Real tp = m_horizontal_radius / speed;
            const amrex::Real ti =
                m_integral_factor * tp * (m_averaging_time + tp);
            inv_ti_l = 1.0_rt / ti;
            inv_tau_l = m_integral_ratio / ti;
        }
        m_level_tau[l] = 1.0_rt / inv_tau_l;
        m_level_ti[l] = (inv_ti_l > 0.0_rt) ? 1.0_rt / inv_ti_l : 0.0_rt;

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            const int idx = (l * AMREX_SPACEDIM) + n;
            if (m_monitor_only ||
                ((n == AMREX_SPACEDIM - 1) && !m_force_vertical)) {
                m_force[idx] = 0.0_rt;
                continue;
            }
            const amrex::Real err = m_level_vel[idx] - m_footprint_mean[idx];
            if (m_gain_schedule) {
                // The gains change in time, so integrate the force itself
                if (advance) {
                    m_integral[idx] += err * dt * inv_ti_l;
                }
                m_force[idx] = (err * inv_tau_l) + m_integral[idx];
            } else {
                if (advance) {
                    m_integral[idx] += err * dt;
                }
                m_force[idx] = (err * inv_tau_l) + (m_integral[idx] * inv_ti_l);
            }
            if ((m_max_force > 0.0_rt) &&
                (std::abs(m_force[idx]) > m_max_force)) {
                // Saturated: cap the force and stop integrating (anti-windup)
                m_force[idx] = std::copysign(m_max_force, m_force[idx]);
                if (advance) {
                    m_integral[idx] -=
                        m_gain_schedule ? err * dt * inv_ti_l : err * dt;
                }
            }
        }
    }
    m_initialized = true;
    copy_to_device(m_force, m_force_d);

    if ((m_output_frequency > 0) &&
        (m_time.time_index() % m_output_frequency == 0)) {
        write_output();
    }
    if (m_time.write_checkpoint()) {
        write_state();
    }
}

void MetMastForcing::write_output() const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    const std::string& dir = m_post_dir;
    const std::string fname = dir + "/metmast_body_force.txt";
    // A fresh run starts a new file, a restart appends to it
    const bool start = !m_output_started &&
                       (m_restart_state.empty() || !amrex::FileExists(fname));
    if (start) {
        amrex::UtilCreateDirectory(dir, 0755);
    }
    std::ofstream ofh(fname, start ? std::ios::out : std::ios::app);
    if (start) {
        ofh << "# step time station level z obs_u obs_v obs_w obs_su obs_sv "
               "obs_sw mean_u mean_v mean_w sigma_u sigma_v sigma_w force_x "
               "force_y force_z tau tau_I\n";
    }
    m_output_started = true;
    ofh << std::setprecision(8);
    for (int s = 0; s < num_stations(); ++s) {
        for (int l = m_level_offset[s]; l < m_level_offset[s + 1]; ++l) {
            ofh << m_time.time_index() << " " << m_time.current_time() << " "
                << s << " " << l - m_level_offset[s] << " " << m_level_z[l];
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                ofh << " " << m_level_vel[(l * AMREX_SPACEDIM) + n];
            }
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                ofh << " " << m_level_sigma[(l * AMREX_SPACEDIM) + n];
            }
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                ofh << " " << footprint_velocity(l, n);
            }
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                ofh << " " << footprint_sigma(l, n);
            }
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                ofh << " " << body_force(l, n);
            }
            ofh << " " << m_level_tau[l] << " " << m_level_ti[l] << "\n";
        }
    }
}

void MetMastForcing::write_state() const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    const std::string& dir = m_post_dir;
    amrex::UtilCreateDirectory(dir, 0755);
    const std::string fname =
        amrex::Concatenate(dir + "/metmast_state", m_time.time_index()) +
        ".txt";
    std::ofstream ofh(fname);
    ofh << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    ofh << num_levels() << "\n";
    for (int idx = 0; idx < num_levels() * AMREX_SPACEDIM; ++idx) {
        ofh << m_footprint_mean[idx] << " " << m_second_moment[idx] << " "
            << m_integral[idx] << "\n";
    }
}

void MetMastForcing::read_state(const std::string& fname) const
{
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("MetMastForcing: cannot find state file " + fname);
    }
    int nlevels = 0;
    if (!(ifh >> nlevels) || (nlevels != num_levels())) {
        amrex::Abort(
            "MetMastForcing: state file " + fname +
            " does not match the station levels");
    }
    for (int idx = 0; idx < num_levels() * AMREX_SPACEDIM; ++idx) {
        if (!(ifh >> m_footprint_mean[idx] >> m_second_moment[idx] >>
              m_integral[idx])) {
            amrex::Abort("MetMastForcing: cannot parse state file " + fname);
        }
    }
    m_initialized = true;
}

void MetMastForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const metmast::StationView st{
        .x = m_station_x_d.data(),
        .y = m_station_y_d.data(),
        .offset = m_level_offset_d.data(),
        .z = m_level_z_d.data(),
        .vel = m_level_vel_d.data(),
        .sigma = m_level_sigma_d.data(),
        .nstations = num_stations()};

    const amrex::Real dt = m_time.delta_t();
    const metmast::ForcingParams prm{
        .inv_rh2 = 1.0_rt / (m_horizontal_radius * m_horizontal_radius),
        .inv_rh2_avg = 1.0_rt / (m_averaging_radius * m_averaging_radius),
        .inv_rz2 = 1.0_rt / (m_vertical_radius * m_vertical_radius),
        .cutoff = m_damping_radius,
        .sigma_factor = m_sigma_factor,
        .inv_tau = 1.0_rt / m_timescale,
        .max_rate = (dt > 0.0_rt) ? 1.0_rt / dt
                                  : std::numeric_limits<amrex::Real>::max(),
        .half_gate = 0.5_rt * m_gate_length};

    const auto& geom = m_mesh.Geom(lev);
    const auto& repo = m_sim.repo();

    if (m_body_force_mode) {
        if (m_time.time_index() != m_last_update_step) {
            update_body_force(st, prm);
            m_last_update_step = m_time.time_index();
        }
        if (m_monitor_only) {
            return;
        }
        if (repo.field_exists("terrain_height")) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                repo.int_field_exists("terrain_blank"),
                "MetMastForcing: terrain_height requires terrain_blank");
            metmast::add_body_force<true>(
                geom, st, prm, m_force_d.data(),
                repo.get_field("terrain_height")(lev).const_arrays(),
                repo.get_int_field("terrain_blank")(lev).const_arrays(),
                src_term);
        } else {
            metmast::add_body_force<false>(
                geom, st, prm, m_force_d.data(),
                amrex::MultiArray4<amrex::Real const>(),
                amrex::MultiArray4<int const>(), src_term);
        }
        return;
    }

    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();

    if (repo.field_exists("terrain_height")) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            repo.int_field_exists("terrain_blank"),
            "MetMastForcing: terrain_height requires terrain_blank");
        metmast::add_forcing<true>(
            geom, st, prm, vel_arrs,
            repo.get_field("terrain_height")(lev).const_arrays(),
            repo.get_int_field("terrain_blank")(lev).const_arrays(), src_term);
    } else {
        metmast::add_forcing<false>(
            geom, st, prm, vel_arrs, amrex::MultiArray4<amrex::Real const>(),
            amrex::MultiArray4<int const>(), src_term);
    }
}

} // namespace kynema_sgf::pde::icns
