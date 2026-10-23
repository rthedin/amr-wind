#ifndef SURFACESTABILITY_H
#define SURFACESTABILITY_H

#include <cmath>
#include <cstdint>
#include <string>
#include "src/utilities/constants.H"
#include "src/utilities/trig_ops.H"
#include "AMReX_REAL.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_Algorithm.H"

using namespace amrex::literals;

namespace kynema_sgf {
class CFDSim;
class Field;
class IntField;
} // namespace kynema_sgf

/** Surface-layer stability shared by the terrain drag source terms
 *  \ingroup we_abl
 *
 *  The stability of the terrain surface layer is set by one of four inputs,
 *  as in ABLWallFunction and the amrTerrain 1-D RANS solver:
 *  the Monin-Obukhov length, the surface heat flux, a surface heating or
 *  cooling rate, or the surface temperature. The functions below return the
 *  Monin-Obukhov length and the surface fluxes at one cell from the wind speed
 *  and potential temperature at a reference height, so that DragForcing,
 *  DragTempForcing and Kosovic use the same surface layer.
 *
 *  The stability functions follow Dyer (1974), as in MOData.
 */
namespace kynema_sgf::surface_stability {

//! Input that sets the surface-layer stability (ABL.surface_stability_type)
enum class StabilityType : std::uint8_t {
    MONIN_OBUKHOV_LENGTH =
        0,        ///< Obukhov length given (ABL.monin_obukhov_length)
    HEAT_FLUX,    ///< Kinematic heat flux given (ABL.surface_temp_flux)
    HEATING_RATE, ///< Surface heating/cooling rate (ABL.surface_temp_rate)
    SURFACE_TEMPERATURE, ///< Surface temperature given (ABL.surface_temp_init)
};

//! Parsed surface-stability inputs, captured by value in device kernels
struct StabilityParams
{
    StabilityType type{StabilityType::MONIN_OBUKHOV_LENGTH};
    //! True when ABL.wall_het_model = "mol"; otherwise the surface is neutral
    bool non_neutral{false};
    amrex::Real obukhov_len{constants::LARGE_NUM}; ///< Obukhov length (m)
    amrex::Real surf_temp_flux{0.0_rt};   ///< Kinematic heat flux (K m/s)
    amrex::Real surf_temp{300.0_rt};      ///< Current surface temperature (K)
    amrex::Real surf_temp_init{300.0_rt}; ///< Surface temperature at tstart (K)
    amrex::Real surf_temp_rate{0.0_rt};   ///< Heating/cooling rate (K/hr)
    amrex::Real surf_temp_rate_tstart{0.0_rt}; ///< Start of heating (s)
    amrex::Real kappa{0.41_rt};
    amrex::Real gravity{9.81_rt};
    amrex::Real gamma_m{5.0_rt};
    amrex::Real beta_m{16.0_rt};
    amrex::Real gamma_h{5.0_rt};
    amrex::Real beta_h{16.0_rt};

    //! True when the Obukhov length depends on the local flow
    [[nodiscard]] AMREX_GPU_HOST_DEVICE bool is_local() const
    {
        return non_neutral && (type != StabilityType::MONIN_OBUKHOV_LENGTH);
    }
};

//! Surface-layer state at one cell
struct SurfaceLayerState
{
    amrex::Real utau{0.0_rt};                      ///< Friction velocity (m/s)
    amrex::Real obukhov_len{constants::LARGE_NUM}; ///< Obukhov length (m)
    amrex::Real surf_temp_flux{0.0_rt};            ///< Heat flux (K m/s)
    amrex::Real surf_temp{0.0_rt};                 ///< Surface temperature (K)
};

//! Momentum stability correction psi_m(zeta)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real psi_m(
    const amrex::Real zeta, const amrex::Real beta_m, const amrex::Real gamma_m)
{
    if (zeta > 0) {
        return -gamma_m * zeta;
    }
    const amrex::Real x = std::sqrt(std::sqrt(1.0_rt - (beta_m * zeta)));
    return (2.0_rt * std::log(0.5_rt * (1.0_rt + x))) +
           std::log(0.5_rt * (1.0_rt + x * x)) - (2.0_rt * std::atan(x)) +
           utils::half_pi();
}

//! Heat stability correction psi_h(zeta)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real psi_h(
    const amrex::Real zeta, const amrex::Real beta_h, const amrex::Real gamma_h)
{
    if (zeta > 0) {
        return -gamma_h * zeta;
    }
    const amrex::Real x = std::sqrt(1 - (beta_h * zeta));
    return 2.0_rt * std::log(0.5_rt * (1 + x));
}

//! Smallest zeta = z/L of the per-cell solve (unstable side)
static constexpr amrex::Real zeta_min = -5.0_rt;
//! Largest zeta = z/L of the per-cell solve (stable side)
static constexpr amrex::Real zeta_max = 2.0_rt;

/** Solve the surface layer at one cell
 *
 *  \param p Stability inputs
 *  \param wspd Horizontal wind speed at the reference height (m/s)
 *  \param theta Potential temperature at the reference height (K)
 *  \param zref Reference height above the surface (m)
 *  \param z0 Roughness length, used for momentum and heat (m); limited to
 *  zref/e
 *
 *  For the Obukhov-length input the result is closed form. For the other
 *  inputs zeta = zref/L is the root of
 *
 *    heat flux:            zeta + zref kappa g q (ln(zref/z0) - psi_m)^3
 *                                 / ((kappa U)^3 theta) = 0
 *    surface temperature:  zeta - Ri_b (ln(zref/z0) - psi_m)^2
 *                                 / (ln(zref/z0) - psi_h) = 0,
 *                          Ri_b = g zref (theta - T_s) / (theta U^2)
 *
 *  found by bisection on [zeta_lo, zeta_hi]. zeta_hi = zeta_max. zeta_lo is
 *  zeta_min or, if larger, the zeta at which ln(zref/z0) - psi_h falls to half
 *  of ln(zref/z0), which keeps both log-law denominators positive. A stable
 *  heat flux has two roots below the critical stability; the search stops at
 *  the first one, the physical root with the larger friction velocity.
 *
 *  When there is no root in the range (a stable surface below the critical
 *  wind speed, or a very unstable one at low wind) zeta is held at the bound.
 *  For the heat-flux input the flux is then reduced to the one the bound can
 *  carry, q = -zeta u_*^3 theta / (zref kappa g), so that the surface
 *  temperature stays bounded; for the surface-temperature inputs the heat
 *  flux follows from the surface temperature.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE SurfaceLayerState solve(
    const StabilityParams& p,
    const amrex::Real wspd,
    const amrex::Real theta,
    const amrex::Real zref,
    const amrex::Real z0)
{
    constexpr int n_bisect = 60;
    // Guard: a roughness length at or above the reference height (for example
    // a forest roughness on a fine mesh) would make ln(zref/z0) zero or
    // negative and the fluxes infinite or of the wrong sign. The roughness is
    // limited to zref/e so that ln(zref/z0) >= 1.
    const amrex::Real log_ratio = std::log(
        zref / amrex::min<amrex::Real>(z0, zref * 0.36787944117144233_rt));
    // Guard: a stagnant cell has zero wind speed, which the heat-flux root
    // divides by
    const amrex::Real speed = amrex::max<amrex::Real>(wspd, 1.0e-6_rt);
    SurfaceLayerState s;

    if (!p.non_neutral) {
        s.utau = p.kappa * wspd / log_ratio;
        s.surf_temp = theta;
        return s;
    }

    if (p.type == StabilityType::MONIN_OBUKHOV_LENGTH) {
        const amrex::Real zeta = zref / p.obukhov_len;
        s.obukhov_len = p.obukhov_len;
        s.utau =
            p.kappa * wspd / (log_ratio - psi_m(zeta, p.beta_m, p.gamma_m));
        const amrex::Real thetastar =
            theta * s.utau * s.utau / (p.kappa * p.gravity * p.obukhov_len);
        s.surf_temp_flux = -s.utau * thetastar;
        s.surf_temp = theta - (thetastar / p.kappa *
                               (log_ratio - psi_h(zeta, p.beta_h, p.gamma_h)));
        return s;
    }

    const bool is_flux = (p.type == StabilityType::HEAT_FLUX);
    // Sign of the stability: q < 0 or T_s < theta is stable (zeta > 0)
    const amrex::Real forcing =
        is_flux ? -p.surf_temp_flux : (theta - p.surf_temp);
    amrex::Real zeta = 0.0_rt;
    bool limited = false;
    if (std::abs(forcing) > constants::EPS) {
        // Residual whose root is zeta; positive below the root on the
        // unstable side and negative below the first root on the stable side
        const amrex::Real ku = p.kappa * speed;
        const amrex::Real c_flux = zref * p.kappa * p.gravity *
                                   p.surf_temp_flux / (ku * ku * ku * theta);
        const amrex::Real rib =
            p.gravity * zref * (theta - p.surf_temp) / (theta * speed * speed);
        auto residual = [=](const amrex::Real z) {
            const amrex::Real dm = log_ratio - psi_m(z, p.beta_m, p.gamma_m);
            if (is_flux) {
                return z + (c_flux * dm * dm * dm);
            }
            const amrex::Real dh = log_ratio - psi_h(z, p.beta_h, p.gamma_h);
            return z - (rib * dm * dm / dh);
        };

        amrex::Real lo = 0.0_rt;
        amrex::Real hi = 0.0_rt;
        if (forcing > 0.0_rt) {
            // Stable. For the heat flux the residual rises to a maximum and
            // falls again; search only up to that maximum so that the first
            // root is found.
            hi = zeta_max;
            if (is_flux && p.gamma_m > 0.0_rt) {
                const amrex::Real z_peak =
                    (1.0_rt / std::sqrt(-3.0_rt * c_flux * p.gamma_m) -
                     log_ratio) /
                    p.gamma_m;
                hi = amrex::min<amrex::Real>(hi, z_peak);
            }
            if ((hi <= 0.0_rt) || (residual(hi) < 0.0_rt)) {
                zeta = zeta_max;
                limited = true;
            } else {
                lo = 0.0_rt;
            }
        } else {
            // Unstable
            const amrex::Real x =
                (2.0_rt * std::exp(0.25_rt * log_ratio)) - 1.0_rt;
            const amrex::Real zeta_h = (1.0_rt - (x * x)) / p.beta_h;
            lo = amrex::max<amrex::Real>(zeta_min, zeta_h);
            hi = 0.0_rt;
            if (residual(lo) > 0.0_rt) {
                zeta = lo;
                limited = true;
            }
        }
        if (!limited) {
            // residual(lo) and residual(hi) have opposite signs
            const amrex::Real r_lo = residual(lo);
            for (int n = 0; n < n_bisect; ++n) {
                const amrex::Real mid = 0.5_rt * (lo + hi);
                if ((residual(mid) > 0.0_rt) == (r_lo > 0.0_rt)) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            zeta = 0.5_rt * (lo + hi);
        }
    }

    const amrex::Real psim = psi_m(zeta, p.beta_m, p.gamma_m);
    const amrex::Real psih = psi_h(zeta, p.beta_h, p.gamma_h);
    s.utau = p.kappa * speed / (log_ratio - psim);
    s.obukhov_len = (zeta != 0.0_rt) ? zref / zeta : constants::LARGE_NUM;
    if (is_flux) {
        s.surf_temp_flux = limited ? -zeta * s.utau * s.utau * s.utau * theta /
                                         (zref * p.kappa * p.gravity)
                                   : p.surf_temp_flux;
        s.surf_temp = theta + (s.surf_temp_flux * (log_ratio - psih) /
                               (s.utau * p.kappa));
    } else {
        s.surf_temp = p.surf_temp;
        s.surf_temp_flux =
            -(theta - s.surf_temp) * s.utau * p.kappa / (log_ratio - psih);
    }
    return s;
}

/** Potential temperature at height z from the surface-layer state
 *
 *  theta(z) = T_s - q (ln(z/z0) - psi_h(z/L)) / (kappa u_*)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real theta_at(
    const StabilityParams& p,
    const SurfaceLayerState& s,
    const amrex::Real z,
    const amrex::Real z0)
{
    const amrex::Real psih =
        p.non_neutral ? psi_h(z / s.obukhov_len, p.beta_h, p.gamma_h) : 0.0_rt;
    return s.surf_temp -
           (s.surf_temp_flux * (std::log(z / z0) - psih) / (p.kappa * s.utau));
}

/** Read the surface-stability inputs from the ABL namespace
 *
 *  \param caller Name of the calling class, used in the abort messages
 *  \param reference_temperature Default for ABL.surface_temp_init with the
 *  heating-rate input
 */
StabilityParams
parse(const std::string& caller, amrex::Real reference_temperature);

/** Update the time-dependent surface temperature of the heating-rate input
 *
 *  T_s(t) = T_s0 + rate max(t - t_start, 0) / 3600
 */
void update(StabilityParams& p, amrex::Real time);

//! Name of a stability type as it appears in the input file
std::string type_name(StabilityType type);

//! Abort unless the stability type is the Obukhov length
void require_obukhov_length(
    const StabilityParams& p, const std::string& caller);

/** Mean of 1/L over the terrain-drag cells of one level
 *
 *  Each terrain-drag cell is solved with the wind and temperature of the cell
 *  above it at 1.5 dz. The inverse Obukhov length is averaged, rather than L
 *  itself, so that cells near neutral do not dominate the mean. Returns zero
 *  when the level has no terrain-drag cells.
 */
amrex::Real mean_inverse_obukhov_length(
    const StabilityParams& p,
    const CFDSim& sim,
    int lev,
    const Field& velocity,
    const Field& temperature);

} // namespace kynema_sgf::surface_stability

#endif /* SURFACESTABILITY_H */
