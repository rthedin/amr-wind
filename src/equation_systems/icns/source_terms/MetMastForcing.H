#ifndef METMASTFORCING_H
#define METMASTFORCING_H

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"

namespace kynema_sgf::pde::icns {

namespace metmast {
struct StationView;
struct ForcingParams;
} // namespace metmast

/** Relax the velocity towards met-mast and lidar measurements
 *  \ingroup icns_src
 *
 *  Each measurement station is a horizontal location and a profile of one or
 *  more heights above the terrain. A met-mast point is a one-level profile and
 *  a lidar is a multi-level profile with the standard deviation of the wind.
 *  Station \f$i\f$ has the weight
 *
 *  \f[ w_i = \exp\left(-\frac{1}{4}\left[\frac{r_i^2}{R_h^2} +
 *  \frac{d_i^2}{R_z^2}\right]\right) \f]
 *
 *  with \f$r_i\f$ the horizontal distance to the station and \f$d_i\f$ the
 *  vertical distance outside the measured height range (zero inside it). The
 *  target is the measured velocity \f$\bar{u}_i\f$ widened to a band of
 *  \f$\pm\alpha \sigma_i\f$, so that only the part of the velocity outside
 *  the measured variability is forced:
 *
 *  \f[ t_i = \min(\max(u, \bar{u}_i - \alpha\sigma_i), \bar{u}_i +
 *  \alpha\sigma_i) \f]
 *
 *  The stations are combined into one relaxation
 *
 *  \f[ S = -\min\left(\frac{W}{\tau}, \frac{1}{\Delta t}\right) (u - \bar{t}),
 *  \quad W = \min\left(1, \sum_i w_i\right), \quad \bar{t} = \frac{\sum_i w_i
 *  t_i}{\sum_i w_i} \f]
 *
 *  The heights are above the local terrain when TerrainDrag is active, and no
 *  forcing is applied in the cells inside the terrain.
 *
 *  The relaxation above acts on the instantaneous velocity, so in LES it also
 *  damps the resolved turbulence. With `ABL.metmast_forcing_type = body_force`,
 * the forcing instead acts on the model's footprint average, as a virtual
 * lidar. Each step, the velocity is averaged over the weights \f$w_i\f$ of
 * every level \f$l\f$ of every station, with the horizontal radius
 *  `ABL.metmast_averaging_radius` (default \f$R_h\f$) to match the lidar
 *  measurement volume, and filtered in time over \f$T_{avg}\f$:
 *
 *  \f[ \bar{U}_{il} \leftarrow \bar{U}_{il} + \frac{\Delta t}{T_{avg}}
 *  \left(\frac{\sum w_{il} u \Delta V}{\sum w_{il} \Delta V} -
 *  \bar{U}_{il}\right) \f]
 *
 *  A proportional-integral controller turns the error \f$e_{il} =
 *  \bar{u}_{il} - \bar{U}_{il}\f$ into a body force
 *
 *  \f[ F_{il} = \frac{e_{il}}{\tau} + \frac{1}{\tau_I} \int e_{il}\,dt \f]
 *
 *  which is interpolated to the cells with the same weights:
 *  \f$S = \sum_i w_i F_i / \max(1, \sum_i w_i)\f$. The force does not depend
 *  on the local velocity, so it shifts the mean without damping the
 *  fluctuations. The vertical component is only forced with
 *  `ABL.metmast_force_vertical = true`.
 *
 *  With `ABL.metmast_gate_length`, each gate averages only the cells within
 *  half the gate length of its height, like a lidar range gate, instead of the
 *  forcing weights and their taper beyond the end gates.
 *  `ABL.metmast_forcing_type = monitor` computes and writes the footprint
 *  averages without applying any force.
 */
class MetMastForcing : public MomentumSource::Register<MetMastForcing>
{
public:
    static std::string identifier() { return "MetMastForcing"; }

    explicit MetMastForcing(const CFDSim& sim);

    ~MetMastForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

    //! Number of measurement stations (met-mast points and lidars)
    [[nodiscard]] int num_stations() const
    {
        return static_cast<int>(m_station_x.size());
    }

    //! Total number of levels over all stations
    [[nodiscard]] int num_levels() const
    {
        return static_cast<int>(m_level_z.size());
    }

    //! Body-force mode filtered footprint velocity of a level (global index)
    [[nodiscard]] amrex::Real footprint_velocity(int ilev, int n) const
    {
        return m_footprint_mean[(ilev * AMREX_SPACEDIM) + n];
    }

    //! Body-force mode footprint standard deviation of a level (global index)
    [[nodiscard]] amrex::Real footprint_sigma(int ilev, int n) const;

    //! Body-force mode body force of a level (global index)
    [[nodiscard]] amrex::Real body_force(int ilev, int n) const
    {
        return m_force[(ilev * AMREX_SPACEDIM) + n];
    }

private:
    //! Average the velocity over the footprints and update the controller
    void update_body_force(
        const metmast::StationView& st,
        const metmast::ForcingParams& prm) const;

    //! Write the footprint averages and forces to the output file
    void write_output() const;

    //! Write the controller state for a restart
    void write_state() const;

    //! Read the controller state written by write_state
    void read_state(const std::string& fname) const;

    //! Read met-mast points, one per line: x y z u v w T
    void read_point_file(const std::string& fname);

    //! Read a lidar profile: x y on the first line, then z u v w su sv sw
    void read_profile_file(const std::string& fname);

    void add_level(
        amrex::Real z,
        const amrex::Array<amrex::Real, AMREX_SPACEDIM>& vel,
        const amrex::Array<amrex::Real, AMREX_SPACEDIM>& sigma);

    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const CFDSim& m_sim;

    //! Station locations
    amrex::Vector<amrex::Real> m_station_x;
    amrex::Vector<amrex::Real> m_station_y;
    //! Levels of station i are [m_level_offset[i], m_level_offset[i+1])
    amrex::Vector<int> m_level_offset{0};
    //! Height above terrain of each level
    amrex::Vector<amrex::Real> m_level_z;
    //! Velocity and its standard deviation, AMREX_SPACEDIM values per level
    amrex::Vector<amrex::Real> m_level_vel;
    amrex::Vector<amrex::Real> m_level_sigma;

    amrex::Gpu::DeviceVector<amrex::Real> m_station_x_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_station_y_d;
    amrex::Gpu::DeviceVector<int> m_level_offset_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_level_z_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_level_vel_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_level_sigma_d;

    //! Relaxation time scale
    amrex::Real m_timescale{30};
    //! Horizontal radius of the forcing footprint
    amrex::Real m_horizontal_radius{500};
    //! Vertical radius of the forcing footprint
    amrex::Real m_vertical_radius{25};
    //! Cutoff on the normalized squared distance
    amrex::Real m_damping_radius{1400};
    //! Half-width of the tolerance band in standard deviations
    amrex::Real m_sigma_factor{1};

    //! Body-force mode (true) or relaxation of the local velocity (false)
    bool m_body_force_mode{false};
    //! Time filter of the footprint averages
    amrex::Real m_averaging_time{120};
    //! Horizontal radius of the footprint average (default: forcing radius)
    amrex::Real m_averaging_radius{500};
    //! Range-gate length of the footprint average (0: forcing weights)
    amrex::Real m_gate_length{0};
    //! Compute and write the footprint averages without forcing
    bool m_monitor_only{false};
    //! Integral time scale of the controller (0 or negative: no integral)
    amrex::Real m_integral_timescale{-1};
    //! Also force the vertical velocity in body-force mode
    bool m_force_vertical{false};
    //! Largest body force per component in body-force mode (0 or negative:
    //! none)
    amrex::Real m_max_force{0};
    //! Steps between outputs of the footprint averages (0: no output)
    int m_output_frequency{10};
    //! Controller state to read on the first update of a restart
    std::string m_restart_state;
    //! Output directory, the same as the sampling output
    std::string m_post_dir{"post_processing"};

    //! Body-force mode state, AMREX_SPACEDIM values per level
    mutable amrex::Vector<amrex::Real> m_footprint_mean;
    mutable amrex::Vector<amrex::Real> m_second_moment;
    mutable amrex::Vector<amrex::Real> m_integral;
    mutable amrex::Vector<amrex::Real> m_force;
    mutable amrex::Gpu::DeviceVector<amrex::Real> m_force_d;
    mutable amrex::Gpu::DeviceVector<amrex::Real> m_sums_d;
    mutable bool m_initialized{false};
    mutable int m_last_update_step{-1};
    mutable bool m_output_started{false};
};

} // namespace kynema_sgf::pde::icns
#endif
