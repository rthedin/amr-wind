#include "src/physics/ForestCanopy.H"
#include "AMReX_ParmParse.H"
#include <cmath>

namespace kynema_sgf::forestdrag {

CanopyTurbulence parse_canopy_turbulence()
{
    CanopyTurbulence canopy;
    amrex::ParmParse pp("ForestDrag");
    pp.query("canopy_tke", canopy.m_enabled);
    pp.query("canopy_beta_p", canopy.m_beta_p);
    pp.query("canopy_beta_d", canopy.m_beta_d);
    if (canopy.m_beta_p < 0.0_rt || canopy.m_beta_d < 0.0_rt) {
        amrex::Abort(
            "ForestDrag: canopy_beta_p and canopy_beta_d must not be "
            "negative");
    }
    return canopy;
}

void add_canopy_tke_source(
    amrex::MultiFab& src_term,
    const amrex::MultiFab& forest_drag,
    const amrex::MultiFab& velocity,
    const amrex::MultiFab& tke,
    const amrex::Real beta_p,
    const amrex::Real beta_d,
    const amrex::Real dt)
{
    auto const& src_arrs = src_term.arrays();
    auto const& forest_arrs = forest_drag.const_arrays();
    auto const& vel_arrs = velocity.const_arrays();
    auto const& tke_arrs = tke.const_arrays();
    amrex::ParallelFor(
        src_term, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real fd = forest_arrs[nbx](i, j, k);
            if (fd > 0.0_rt) {
                const amrex::Real ux = vel_arrs[nbx](i, j, k, 0);
                const amrex::Real uy = vel_arrs[nbx](i, j, k, 1);
                const amrex::Real uz = vel_arrs[nbx](i, j, k, 2);
                const amrex::Real ws =
                    std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
                // expm1 is in (-1, 0]: the sink never exceeds k / dt
                src_arrs[nbx](i, j, k) +=
                    (beta_p * fd * ws * ws * ws) +
                    (tke_arrs[nbx](i, j, k) *
                     std::expm1(-beta_d * fd * ws * dt) / dt);
            }
        });
}

} // namespace kynema_sgf::forestdrag
