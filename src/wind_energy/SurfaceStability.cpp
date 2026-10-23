#include <cmath>
#include "src/wind_energy/SurfaceStability.H"
#include "src/CFDSim.H"
#include "src/core/Field.H"
#include "src/core/IntField.H"
#include "src/utilities/trig_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::surface_stability {

namespace {
void get_required(
    const amrex::ParmParse& pp,
    const std::string& caller,
    const std::string& key,
    const std::string& type,
    amrex::Real& value)
{
    if (!pp.contains(key)) {
        amrex::Abort(
            caller + ": ABL." + key + " is required with " +
            "ABL.surface_stability_type = " + type);
    }
    pp.get(key, value);
}

void check_finite(
    const std::string& caller, const std::string& key, amrex::Real value)
{
    if (!std::isfinite(value)) {
        amrex::Abort(caller + ": ABL." + key + " must be finite");
    }
}
} // namespace

std::string type_name(const StabilityType type)
{
    switch (type) {
    case StabilityType::MONIN_OBUKHOV_LENGTH:
        return "monin_obukhov_length";
    case StabilityType::HEAT_FLUX:
        return "heat_flux";
    case StabilityType::HEATING_RATE:
        return "heating_rate";
    case StabilityType::SURFACE_TEMPERATURE:
        return "surface_temperature";
    }
    return "unknown";
}

StabilityParams
parse(const std::string& caller, const amrex::Real reference_temperature)
{
    StabilityParams p;
    amrex::ParmParse pp("ABL");

    std::string wall_het_model{"none"};
    pp.query("wall_het_model", wall_het_model);
    p.non_neutral = (wall_het_model == "mol");

    std::string type{"monin_obukhov_length"};
    pp.query("surface_stability_type", type);
    if (type == "monin_obukhov_length") {
        p.type = StabilityType::MONIN_OBUKHOV_LENGTH;
    } else if (type == "heat_flux") {
        p.type = StabilityType::HEAT_FLUX;
    } else if (type == "heating_rate") {
        p.type = StabilityType::HEATING_RATE;
    } else if (type == "surface_temperature") {
        p.type = StabilityType::SURFACE_TEMPERATURE;
    } else {
        amrex::Abort(
            caller + ": ABL.surface_stability_type = " + type +
            " is not one of monin_obukhov_length, heat_flux, heating_rate, "
            "surface_temperature");
    }
    if ((p.type != StabilityType::MONIN_OBUKHOV_LENGTH) && !p.non_neutral) {
        amrex::Abort(
            caller + ": ABL.surface_stability_type = " + type +
            " requires ABL.wall_het_model = mol");
    }

    pp.query("kappa", p.kappa);
    pp.query("mo_gamma_m", p.gamma_m);
    pp.query("mo_beta_m", p.beta_m);
    pp.query("mo_gamma_h", p.gamma_h);
    pp.query("mo_beta_h", p.beta_h);
    if (p.kappa <= 0.0_rt) {
        amrex::Abort(caller + ": ABL.kappa must be > 0");
    }
    if ((p.gamma_m < 0.0_rt) || (p.beta_m < 0.0_rt) || (p.gamma_h < 0.0_rt) ||
        (p.beta_h < 0.0_rt)) {
        amrex::Abort(
            caller +
            ": ABL.mo_gamma_m, ABL.mo_beta_m, ABL.mo_gamma_h and "
            "ABL.mo_beta_h must be >= 0");
    }

    amrex::Vector<amrex::Real> gravity{0.0_rt, 0.0_rt, -9.81_rt};
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", gravity);
    p.gravity = std::abs(gravity[2]);

    // The flat-wall ABLWallFunction reads the same surface temperature keys
    // and picks one of them in a fixed order. Any key that is not the input of
    // the selected type would make the flat wall and the terrain disagree.
    {
        amrex::Vector<std::string> conflicting;
        switch (p.type) {
        case StabilityType::MONIN_OBUKHOV_LENGTH:
            break;
        case StabilityType::HEAT_FLUX:
            conflicting = {
                "surface_temp_rate", "surface_temp_timetable",
                "surface_temp_flux_timetable", "near_surface_temp_timetable"};
            break;
        case StabilityType::HEATING_RATE:
            conflicting = {
                "surface_temp_flux", "surface_temp_timetable",
                "surface_temp_flux_timetable", "near_surface_temp_timetable"};
            break;
        case StabilityType::SURFACE_TEMPERATURE:
            conflicting = {
                "surface_temp_flux", "surface_temp_rate",
                "surface_temp_timetable", "surface_temp_flux_timetable",
                "near_surface_temp_timetable"};
            break;
        }
        for (const auto& key : conflicting) {
            if (pp.contains(key)) {
                amrex::Abort(
                    caller + ": ABL." + key +
                    " cannot be used with ABL.surface_stability_type = " +
                    type);
            }
        }
    }

    switch (p.type) {
    case StabilityType::MONIN_OBUKHOV_LENGTH:
        pp.query("monin_obukhov_length", p.obukhov_len);
        if (p.non_neutral && (std::abs(p.obukhov_len) <= constants::EPS)) {
            amrex::Abort(
                caller +
                ": ABL.monin_obukhov_length must satisfy |L| > machine "
                "epsilon");
        }
        break;
    case StabilityType::HEAT_FLUX:
        get_required(pp, caller, "surface_temp_flux", type, p.surf_temp_flux);
        check_finite(caller, "surface_temp_flux", p.surf_temp_flux);
        if (std::abs(p.surf_temp_flux) > 10.0_rt) {
            amrex::Abort(
                caller + ": ABL.surface_temp_flux must be in [-10, 10] K m/s");
        }
        break;
    case StabilityType::HEATING_RATE:
        get_required(pp, caller, "surface_temp_rate", type, p.surf_temp_rate);
        check_finite(caller, "surface_temp_rate", p.surf_temp_rate);
        p.surf_temp_init = reference_temperature;
        pp.query("surface_temp_init", p.surf_temp_init);
        pp.query("surface_temp_rate_tstart", p.surf_temp_rate_tstart);
        check_finite(caller, "surface_temp_init", p.surf_temp_init);
        check_finite(
            caller, "surface_temp_rate_tstart", p.surf_temp_rate_tstart);
        if (p.surf_temp_init <= 0.0_rt) {
            amrex::Abort(caller + ": ABL.surface_temp_init must be > 0 K");
        }
        if (p.surf_temp_rate_tstart < 0.0_rt) {
            amrex::Abort(
                caller + ": ABL.surface_temp_rate_tstart must be >= 0");
        }
        break;
    case StabilityType::SURFACE_TEMPERATURE:
        get_required(pp, caller, "surface_temp_init", type, p.surf_temp_init);
        check_finite(caller, "surface_temp_init", p.surf_temp_init);
        if (p.surf_temp_init <= 0.0_rt) {
            amrex::Abort(caller + ": ABL.surface_temp_init must be > 0 K");
        }
        break;
    }
    update(p, 0.0_rt);
    return p;
}

void update(StabilityParams& p, const amrex::Real time)
{
    if (p.type == StabilityType::HEATING_RATE) {
        p.surf_temp =
            p.surf_temp_init +
            (p.surf_temp_rate *
             amrex::max<amrex::Real>(time - p.surf_temp_rate_tstart, 0.0_rt) /
             3600.0_rt);
    } else if (p.type == StabilityType::SURFACE_TEMPERATURE) {
        p.surf_temp = p.surf_temp_init;
    }
}

void require_obukhov_length(const StabilityParams& p, const std::string& caller)
{
    if (p.type != StabilityType::MONIN_OBUKHOV_LENGTH) {
        amrex::Abort(
            caller +
            " uses ABL.monin_obukhov_length for the friction velocity and "
            "ABL.surface_temp_flux for the TKE boundary condition; "
            "ABL.surface_stability_type = " +
            type_name(p.type) + " is not supported");
    }
}

amrex::Real mean_inverse_obukhov_length(
    const StabilityParams& p,
    const CFDSim& sim,
    const int lev,
    const Field& velocity,
    const Field& temperature)
{
    const auto& repo = sim.repo();
    if (!repo.int_field_exists("terrain_drag")) {
        return 0.0_rt;
    }
    const auto& drag = repo.get_int_field("terrain_drag")(lev);
    const auto& drag_arrs = drag.const_arrays();
    const auto& z0_arrs = repo.get_field("terrainz0")(lev).const_arrays();
    const auto& vel_arrs = velocity(lev).const_arrays();
    const auto& theta_arrs = temperature(lev).const_arrays();
    const amrex::Real dz = sim.mesh().Geom(lev).CellSize(2);
    const amrex::Real z0_min = 1.0e-4_rt;

    auto [inv_l_sum, n_cells] = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
        amrex::TypeList<amrex::Real, amrex::Real>{}, drag, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
            -> amrex::GpuTuple<amrex::Real, amrex::Real> {
            if (drag_arrs[nbx](i, j, k) != 1) {
                return {0.0_rt, 0.0_rt};
            }
            const amrex::Real ux = vel_arrs[nbx](i, j, k + 1, 0);
            const amrex::Real uy = vel_arrs[nbx](i, j, k + 1, 1);
            const amrex::Real z0 =
                amrex::max<amrex::Real>(z0_arrs[nbx](i, j, k), z0_min);
            const auto s = solve(
                p, std::sqrt((ux * ux) + (uy * uy)),
                theta_arrs[nbx](i, j, k + 1), 1.5_rt * dz, z0);
            return {1.0_rt / s.obukhov_len, 1.0_rt};
        });
    amrex::Real sums[2] = {inv_l_sum, n_cells};
    amrex::ParallelDescriptor::ReduceRealSum(sums, 2);
    return (sums[1] > 0.0_rt) ? sums[0] / sums[1] : 0.0_rt;
}

} // namespace kynema_sgf::surface_stability
