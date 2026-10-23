#include <limits>
#include "src/wind_energy/MOData.H"
#include "src/wind_energy/SurfaceStability.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

/*
 * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
 * model consistent with Monin-Obukhov similarity theory. Wind Energy,
 * 20(3), 479–489. https://doi.org/10.1002/we.2017
 *
 * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
 * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
 * Demetri P Lalas and Corrado F Ratto, January 1996,
 * https://doi.org/10.1142/2975.
 */

amrex::Real MOData::calc_psi_m(
    const amrex::Real zeta, const amrex::Real beta_m, const amrex::Real gamma_m)
{
    return surface_stability::psi_m(zeta, beta_m, gamma_m);
}

amrex::Real MOData::calc_psi_m(amrex::Real zeta) const
{
    return calc_psi_m(zeta, beta_m, gamma_m);
}

amrex::Real MOData::calc_psi_h(
    const amrex::Real zeta, const amrex::Real beta_h, const amrex::Real gamma_h)
{
    return surface_stability::psi_h(zeta, beta_h, gamma_h);
}

amrex::Real MOData::calc_psi_h(amrex::Real zeta) const
{
    return calc_psi_h(zeta, beta_h, gamma_h);
}

void MOData::update_fluxes(int max_iters)
{
    constexpr amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real zeta = 0.0_rt;
    amrex::Real utau_iter = 0.0_rt;

    // Initialize variables
    amrex::Real psi_m = 0.0_rt;
    amrex::Real psi_h = 0.0_rt;
    utau = kappa * vmag_mean / (std::log(zref / z0));

    int iter = 0;
    do {
        utau_iter = utau;
        switch (alg_type) {
        case ThetaCalcType::HEAT_FLUX:
            surf_temp =
                (surf_temp_flux * (alpha_h * (std::log(zref / z0t) - psi_h)) /
                 (utau * kappa)) +
                theta_mean;
            break;

        case ThetaCalcType::SURFACE_TEMPERATURE:
            surf_temp_flux = -(theta_mean - surf_temp) * utau * kappa /
                             (alpha_h * (std::log(zref / z0t) - psi_h));
            break;

        case ThetaCalcType::NEAR_SURFACE_TEMPERATURE: {
            // theta(z) = theta_s - q phi_h(z) / (kappa utau) holds at both
            // zref and near_surf_height: solve the pair for q and theta_s
            const amrex::Real phi_ref =
                alpha_h * (std::log(zref / z0t) - psi_h);
            const amrex::Real phi_near =
                alpha_h * (std::log(near_surf_height / z0t) -
                           calc_psi_h(zeta * near_surf_height / zref));
            surf_temp_flux = utau * kappa * (near_surf_temp - theta_mean) /
                             (phi_ref - phi_near);
            surf_temp =
                near_surf_temp + (surf_temp_flux * phi_near / (utau * kappa));
            break;
        }
        }

        if (std::abs(surf_temp_flux) > eps) {
            // Stable and unstable ABL conditions
            obukhov_len = -utau * utau * utau * theta_mean /
                          (kappa * gravity * surf_temp_flux);
            zeta = zref / obukhov_len;
        } else {
            // Neutral conditions
            obukhov_len = std::numeric_limits<amrex::Real>::max();
            zeta = 0.0_rt;
        }
        psi_m = calc_psi_m(zeta);
        psi_h = calc_psi_h(zeta);
        utau = kappa * vmag_mean / (std::log(zref / z0) - psi_m);
        ++iter;
    } while ((std::abs(utau_iter - utau) > 1.0e-5_rt) && iter <= max_iters);

    if (iter >= max_iters) {
        amrex::Print()
            << "MOData::update_fluxes: Convergence criteria not met after "
            << max_iters << " iterations\nObuhov length = " << obukhov_len
            << " zeta = " << zeta << "\npsi_m = " << psi_m
            << " psi_h = " << psi_h << "\nutau = " << utau
            << " Tsurf = " << surf_temp << " q = " << surf_temp_flux << '\n';
    }
}

} // namespace kynema_sgf
