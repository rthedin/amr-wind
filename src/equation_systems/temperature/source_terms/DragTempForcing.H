#ifndef DRAGTEMPFORCING_H
#define DRAGTEMPFORCING_H

#include "src/equation_systems/temperature/TemperatureSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "src/utilities/constants.H"
#include "src/wind_energy/SurfaceStability.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

class DragTempForcing : public TemperatureSource::Register<DragTempForcing>
{
public:
    static std::string identifier() { return "DragTempForcing"; }

    explicit DragTempForcing(const CFDSim& sim);

    ~DragTempForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field& m_temperature;
    amrex::Real m_drag_coefficient{1.0_rt};

    //! Surface-layer stability inputs (ABL.surface_stability_type)
    surface_stability::StabilityParams m_stability;
    amrex::Real m_soil_temperature{300};
    amrex::Real m_forcing_time_factor{5.0_rt};
};

} // namespace kynema_sgf::pde::temperature
#endif
