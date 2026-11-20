#ifndef KLAXELL_H
#define KLAXELL_H

#include <cmath>
#include <string>
#include "src/core/ScratchField.H"
#include "src/turbulence/TurbModelBase.H"
#include "src/physics/ForestCanopy.H"
#include "src/wind_energy/SurfaceStability.H"
#include "AMReX_REAL.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_Extension.H"

using namespace amrex::literals;

namespace kynema_sgf::turbulence {

/** Number of cells of a column that TerrainDrag blanks below height h
 *
 *  TerrainDrag blanks cell m when zlo + (m + 1/2) dz <= h. The rounded
 *  estimate floor((h - zlo) / dz + 1/2) can be one off in round-off when h
 *  is at or next to a cell center; it is corrected with the same test, so
 *  that the count matches the blanking.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int terrain_blanked_cells(
    const amrex::Real h, const amrex::Real zlo, const amrex::Real dz) noexcept
{
    int n = static_cast<int>(std::floor(((h - zlo) / dz) + 0.5_rt));
    n = (n > 0) ? n : 0;
    if ((n > 0) && (zlo + ((n - 1 + 0.5_rt) * dz) > h)) {
        --n;
    }
    if (zlo + ((n + 0.5_rt) * dz) <= h) {
        ++n;
    }
    return n;
}

/** Base class for 1-Equation RANS TKE turbulence model
 *  \ingroup turb_model
 */
template <typename Transport>
class KLAxell : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "KLAxell-" + Transport::identifier();
    }

    explicit KLAxell(CFDSim& sim);

    [[nodiscard]] std::string model_name() const override { return "KLAxell"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    //! Do any post advance work
    void post_advance_work() override;

    //! Actions to perform post initialization
    void post_init_actions() override;

    //! Actions to perform post regrid
    void post_regrid_actions() override;

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Update the effective scalar diffusivity field
    void update_scalar_diff(Field& deff, const std::string& name) override;

    /** With TerrainDrag and ``KLAxell.terrain_wall_stencil``: strain rate of
     *  the fluid cells next to the blanked cells from the stencil that flat
     *  ground uses at its bottom wall, the wall being the face of the blanked
     *  cell. Public because it launches a device lambda.
     *
     *  @param velocity Velocity field of the state being updated
     */
    void terrain_wall_strain_rate(const Field& velocity);

    /** With TerrainDrag and ``KLAxell.terrain_blanked_face_length``: height
     *  of the top face of the blanked column, from which the mixing length is
     *  measured in place of the terrain height. Public because it launches a
     *  device lambda.
     *
     *  @param face_height Scratch field (one component) set to the face height
     */
    void terrain_blanked_face_height(ScratchField& face_height) const;

    /** With TerrainDrag and ``KLAxell.terrain_face_stress``: viscosity of the
     *  drag cells sized so that the face above each one carries the wall
     *  stress u*^2 of the DragForcing wall law. Reads mu_{k+1}, so it is
     *  called after the fill of the viscosity. Public because it launches a
     *  device lambda.
     *
     *  @param lev Level to update
     *  @param velocity Velocity field of the state being updated
     *  @param density Density field of the state being updated
     */
    void
    terrain_face_stress(int lev, const Field& velocity, const Field& density);

    /** With TerrainDrag, ``KLAxell.terrain_face_heat_flux`` and
     *  ``ABL.wall_het_model = mol``: heat diffusivity of the drag cells sized
     *  so that the face above each one carries the surface heat flux
     *  -u* theta* given by the Obukhov length. Reads alpha_{k+1}, so it is
     *  called after the fill of alphaeff. Public because it launches a
     *  device lambda.
     *
     *  @param lev Level to update
     *  @param alphaeff Effective thermal diffusivity being updated
     *  @param lam_alpha Laminar thermal diffusivity
     */
    void terrain_face_heat_flux(
        int lev, Field& alphaeff, const ScratchField& lam_alpha);

    //! terrain_face_heat_flux on every level, followed by the fill of
    //! alphaeff
    void terrain_face_heat_flux(Field& alphaeff, const ScratchField& lam_alpha);

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    Field& m_vel;
    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_rho;
    Field* m_tke{nullptr};
    //! Turbulence constant
    amrex::Real m_Cmu{0.556_rt};
    amrex::Real m_Cb_stable{0.25_rt};
    amrex::Real m_Cb_unstable{0.35_rt};
    //! von Karman constant (ABL.kappa), shared with KransAxell
    amrex::Real m_kappa{0.41_rt};
    Field& m_temperature;
    //! Gravity vector (m/s^2)
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
    amrex::Real m_surf_flux{0};
    //! ABL.meso_sponge_start, same default as KransAxell
    amrex::Real m_meso_sponge_start;
    //! KLAxell_coeffs.legacy_alphaeff: recompute Rt in update_alphaeff
    bool m_legacy_alphaeff{false};

    std::unique_ptr<ScratchField> m_gradT;
    //! Pr_t(Rt) from the last update_turbulent_viscosity, used by
    //! update_alphaeff
    std::unique_ptr<ScratchField> m_prandtl_rt;

    //! Canopy length scale limit (ForestDrag.canopy_tke, canopy_length_alpha)
    forestdrag::CanopyTurbulence m_canopy;

    //! Surface-layer inputs (ABL.wall_het_model, ABL.monin_obukhov_length,
    //! ABL.kappa, ABL.mo_beta_m, ABL.mo_gamma_m, incflo.gravity), parsed as
    //! DragForcing and DragTempForcing parse them; used by the face stress
    //! and heat flux below
    surface_stability::StabilityParams m_stability;

    //! TerrainDrag options (``KLAxell.<name>``), all off by default and all
    //! on with ``TerrainDrag.wall_treatment = improved``
    //! Flat-ground wall stencil for the strain rate of the wall cells
    bool m_terrain_wall_stencil{false};
    //! Mixing length measured from the top face of the blanked column
    bool m_terrain_blanked_face_length{false};
    //! Drag-cell viscosity from the wall stress on the face above it
    bool m_terrain_face_stress{false};
    //! Drag-cell heat diffusivity from the surface heat flux on the face
    //! above it (mol only)
    bool m_terrain_face_heat_flux{false};
};

} // namespace kynema_sgf::turbulence

#endif /* KLAXELL_H */
