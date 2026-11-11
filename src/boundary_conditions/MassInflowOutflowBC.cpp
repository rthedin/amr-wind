#include "src/boundary_conditions/MassInflowOutflowBC.H"
#include "AMReX_ParmParse.H"

using namespace amrex::literals;

namespace kynema_sgf {

MassInflowOutflowBC::MassInflowOutflowBC(Field& field, amrex::Orientation ori)
    : m_field(field), m_ori(ori)
{
    // BC bug study 2026-10-07
    amrex::ParmParse pp("incflo");
    pp.query("inout_consistent_ghosts", m_consistent_ghosts);
    amrex::ParmParse pp_eq("ICNS");
    pp_eq.query("clip_inout_backflow", m_clip_backflow);
}

void MassInflowOutflowBC::operator()(
    Field& field, const FieldState /*rho_state*/)
{
    // BC bug study 2026-10-07: act on the field (time state) this op is called
    // on, with the velocity of the same state, into ng normal ghost layers
    const auto& repo = field.repo();
    const auto& velocity =
        repo.get_field("velocity").state(field.field_state());
    const int ncomp = field.num_comp();
    const int idim = m_ori.coordDir();
    const int ng = m_consistent_ghosts ? field.num_grow()[idim] : 1;
    // BC bug study 2026-10-07: clip the normal velocity ghost to no re-entry
    const bool clip = m_clip_backflow && (field.base_name() == "velocity");
    const auto islow = m_ori.isLow();
    const auto ishigh = m_ori.isHigh();
    const int nlevels = m_field.repo().num_active_levels();
    const amrex::IntVect iv_dir = {
        static_cast<int>(idim == 0), static_cast<int>(idim == 1),
        static_cast<int>(idim == 2)};

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& domain = repo.mesh().Geom(lev).Domain();

        amrex::MFItInfo mfi_info{};
        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(field(lev), mfi_info); mfi.isValid(); ++mfi) {
            auto bx = mfi.validbox();
            bx.grow(
                {static_cast<int>(idim != 0), static_cast<int>(idim != 1),
                 static_cast<int>(idim != 2)});
            const auto& bc_a = field(lev).array(mfi);
            const auto& vel = velocity(lev).array(mfi);

            if (islow && (bx.smallEnd(idim) == domain.smallEnd(idim))) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        const amrex::IntVect iv = {i, j, k};
                        const amrex::IntVect ivm = iv - iv_dir;
                        if (vel(ivm[0], ivm[1], ivm[2], idim) < 0) {
                            for (int g = 0; g < ng; g++) {
                                const amrex::IntVect ivg = ivm - g * iv_dir;
                                for (int n = 0; n < ncomp; n++) {
                                    bc_a(ivg[0], ivg[1], ivg[2], n) =
                                        bc_a(i, j, k, n);
                                }
                                if (clip) {
                                    bc_a(ivg[0], ivg[1], ivg[2], idim) =
                                        amrex::min(
                                            bc_a(ivg[0], ivg[1], ivg[2], idim),
                                            0.0_rt);
                                }
                            }
                        }
                    });
            }

            if (ishigh && (bx.bigEnd(idim) == domain.bigEnd(idim))) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        const amrex::IntVect iv = {i, j, k};
                        const amrex::IntVect ivm = iv - iv_dir;
                        if (vel(i, j, k, idim) > 0) {
                            for (int g = 0; g < ng; g++) {
                                const amrex::IntVect ivg = iv + g * iv_dir;
                                for (int n = 0; n < ncomp; n++) {
                                    bc_a(ivg[0], ivg[1], ivg[2], n) =
                                        bc_a(ivm[0], ivm[1], ivm[2], n);
                                }
                                if (clip) {
                                    bc_a(ivg[0], ivg[1], ivg[2], idim) =
                                        amrex::max(
                                            bc_a(ivg[0], ivg[1], ivg[2], idim),
                                            0.0_rt);
                                }
                            }
                        }
                    });
            }
        }
    }
}

} // namespace kynema_sgf
