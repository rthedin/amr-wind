#include "src/equation_systems/temperature/source_terms/DragTempForcing.H"
#include "src/utilities/IOManager.H"
#include "src/transport_models/TransportModel.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

DragTempForcing::DragTempForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("DragTempForcing");
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("soil_temperature", m_soil_temperature);
    pp.query("bc_forcing_time_factor", m_forcing_time_factor);
    m_stability = surface_stability::parse(
        identifier(), sim.transport_model().reference_temperature());
}

DragTempForcing::~DragTempForcing() = default;

void DragTempForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const auto& temperature_arrs =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const bool is_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    if (!is_terrain) {
        amrex::Abort("Need terrain blanking variable to use this source term");
    }
    const auto& blank_arrs =
        this->m_sim.repo().get_int_field("terrain_blank")(lev).const_arrays();
    const auto& drag_arrs =
        this->m_sim.repo().get_int_field("terrain_drag")(lev).const_arrays();
    const auto& terrainz0_arrs =
        this->m_sim.repo().get_field("terrainz0")(lev).const_arrays();
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const amrex::Real drag_coefficient = m_drag_coefficient / dx[2];
    // With the Obukhov-length input the stability corrections are the same in
    // every cell; with the other inputs they are solved per terrain-drag cell
    auto stability = m_stability;
    surface_stability::update(stability, m_time.current_time());
    const bool local_stability = stability.is_local();
    const amrex::Real gravity_mod = stability.gravity;
    const amrex::Real kappa = stability.kappa;
    const amrex::Real z0_min = 1.0e-4_rt;
    const amrex::Real monin_obukhov_length = stability.obukhov_len;
    const auto& dt = m_time.delta_t();
    const amrex::Real psi_m = stability.non_neutral
                                  ? surface_stability::psi_m(
                                        1.5_rt * dx[2] / monin_obukhov_length,
                                        stability.beta_m, stability.gamma_m)
                                  : 0.0_rt;
    const amrex::Real psi_h_neighbour =
        stability.non_neutral ? surface_stability::psi_h(
                                    1.5_rt * dx[2] / monin_obukhov_length,
                                    stability.beta_h, stability.gamma_h)
                              : 0.0_rt;
    const amrex::Real psi_h_cell =
        stability.non_neutral ? surface_stability::psi_h(
                                    0.5_rt * dx[2] / monin_obukhov_length,
                                    stability.beta_h, stability.gamma_h)
                              : 0.0_rt;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real cd_max = 10.0_rt;
    const amrex::Real T0 = m_soil_temperature;
    const amrex::Real time_factor = m_forcing_time_factor;

    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            const auto& vel = vel_arrs[nbx];
            const auto& temperature = temperature_arrs[nbx];
            const auto& blank = blank_arrs[nbx];
            const auto& drag = drag_arrs[nbx];
            const auto& terrainz0 = terrainz0_arrs[nbx];

            const amrex::Real z0 =
                amrex::max<amrex::Real>(terrainz0(i, j, k), z0_min);
            const amrex::Real ux1 = vel(i, j, k, 0);
            const amrex::Real uy1 = vel(i, j, k, 1);
            const amrex::Real uz1 = vel(i, j, k, 2);
            const amrex::Real theta = temperature(i, j, k, 0);
            const amrex::Real theta2 = temperature(i, j, k + 1, 0);
            amrex::Real tTarget = theta;
            if (local_stability && drag(i, j, k) == 1) {
                // Solve the surface layer from the cell above at 1.5 dz and
                // evaluate the profile at the cell centre
                const amrex::Real ux2 = vel(i, j, k + 1, 0);
                const amrex::Real uy2 = vel(i, j, k + 1, 1);
                const auto sl = surface_stability::solve(
                    stability, std::sqrt((ux2 * ux2) + (uy2 * uy2)), theta2,
                    1.5_rt * dx[2], z0);
                tTarget = surface_stability::theta_at(
                    stability, sl, 0.5_rt * dx[2], z0);
            } else if (!local_stability) {
                // The log law is evaluated at 1.5 dz on purpose: the
                // log-law height 0.5 dz of the terrain cell is affected by
                // the terrain forcing
                const amrex::Real wspd = std::sqrt((ux1 * ux1) + (uy1 * uy1));
                const amrex::Real ustar =
                    wspd * kappa / (std::log(1.5_rt * dx[2] / z0) - psi_m);
                const amrex::Real thetastar =
                    theta * ustar * ustar /
                    (kappa * gravity_mod * monin_obukhov_length);
                const amrex::Real surf_temp =
                    theta2 -
                    (thetastar / kappa *
                     (std::log(1.5_rt * dx[2] / z0) - psi_h_neighbour));
                tTarget =
                    surf_temp + (thetastar / kappa *
                                 (std::log(0.5_rt * dx[2] / z0) - psi_h_cell));
            }
            const amrex::Real bc_forcing_t =
                -(tTarget - theta) / (time_factor * dt);
            const amrex::Real m =
                std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
            const amrex::Real Cd = amrex::min<amrex::Real>(
                drag_coefficient / (m + tiny), cd_max / dx[2]);
            src_arrs[nbx](i, j, k, 0) -=
                ((Cd * (theta - T0) * blank(i, j, k, 0)) +
                 (bc_forcing_t * drag(i, j, k)));
        });
}

} // namespace kynema_sgf::pde::temperature
