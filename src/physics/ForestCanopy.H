#ifndef FORESTCANOPY_H
#define FORESTCANOPY_H

#include "AMReX_MultiFab.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::forestdrag {

/** Canopy turbulence inputs of the ForestDrag section
 *
 *  Read by the KransAxell TKE source of the KLAxell model. The
 *  canopy terms act where the forest_drag field (C_d * LAD) is positive, so
 *  they only need ForestDrag in canopy mode.
 */
struct CanopyTurbulence
{
    //! ForestDrag.canopy_tke: add the canopy TKE terms (off by default)
    bool m_enabled{false};
    //! ForestDrag.canopy_beta_p: wake production coefficient
    amrex::Real m_beta_p{1.0_rt};
    //! ForestDrag.canopy_beta_d: short-circuit dissipation coefficient
    amrex::Real m_beta_d{4.0_rt};
};

//! Read the canopy turbulence inputs from the ForestDrag section
CanopyTurbulence parse_canopy_turbulence();

/** Add the canopy TKE source to a TKE source term
 *
 *  S_k = f (beta_p |U|^3 - beta_d |U| k) with f = C_d LAD. The sink is
 *  integrated exactly over the time step, -k (1 - exp(-beta_d f |U| dt)) / dt,
 *  so it removes at most the local k in one step. Only cells with f > 0 are
 *  changed.
 *
 *  \param src_term TKE source term, valid cells updated
 *  \param forest_drag Forest drag field f = C_d LAD (1/m)
 *  \param velocity Velocity field
 *  \param tke Turbulent kinetic energy field
 *  \param beta_p Wake production coefficient
 *  \param beta_d Short-circuit dissipation coefficient
 *  \param dt Time step
 */
void add_canopy_tke_source(
    amrex::MultiFab& src_term,
    const amrex::MultiFab& forest_drag,
    const amrex::MultiFab& velocity,
    const amrex::MultiFab& tke,
    amrex::Real beta_p,
    amrex::Real beta_d,
    amrex::Real dt);

} // namespace kynema_sgf::forestdrag

#endif
