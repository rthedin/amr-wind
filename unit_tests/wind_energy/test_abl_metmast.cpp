#include <cmath>
#include <fstream>
#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/equation_systems/icns/icns.H"
#include "src/equation_systems/icns/source_terms/MetMastForcing.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void write_file(const std::string& fname, const std::string& content)
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << content;
    }
    amrex::ParallelDescriptor::Barrier();
}

// Relaxation expected at a normalized squared distance ri2 from a station
amrex::Real expected_source(
    const amrex::Real ri2,
    const amrex::Real target,
    const amrex::Real vel,
    const amrex::Real tau)
{
    return std::exp(-0.25_rt * ri2) / tau * (target - vel);
}

} // namespace

// Cells are 50 m wide and 25 m tall; centers at x, y = 25 + 50 i and
// z = 12.5 + 25 k
class MetMastTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{16, 16, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("max_grid_size", 8);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{800.0_rt, 800.0_rt, 400.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
    }

    void setup_sim()
    {
        populate_parameters();
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().init_physics();
        sim().repo().get_field("velocity").setVal(m_vel);
    }

    kynema_sgf::Field& src_term()
    {
        return sim().pde_manager().icns().fields().src_term;
    }

    void evaluate(const kynema_sgf::pde::icns::MetMastForcing& forcing)
    {
        src_term().setVal(0.0_rt);
        forcing(0, kynema_sgf::FieldState::New, src_term()(0));
    }

    amrex::Real probe(const int i, const int j, const int k, const int n)
    {
        return utils::field_probe(src_term(), 0, i, j, k, n);
    }

    // The body-force mode averages the old state, the forcing reads neither
    void set_velocity(const amrex::Vector<amrex::Real>& vel)
    {
        auto& velocity = sim().repo().get_field("velocity");
        velocity.setVal(vel);
        velocity.state(kynema_sgf::FieldState::Old).setVal(vel);
    }

    void next_step()
    {
        auto& time = sim().time();
        time.new_timestep();
        time.set_current_cfl(2.0_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
    }

    // One lidar gate at the cell (8, 8, 3) center, in body-force mode
    static void setup_mean_gate(const std::string& fname)
    {
        write_file(fname, "425 425\n87.5 7 2 0 1 1 1\n");
        amrex::ParmParse pp("ABL");
        amrex::Vector<std::string> files{fname};
        pp.addarr("metmast_profile_files", files);
        pp.add("metmast_forcing_type", std::string("body_force"));
        pp.add("metmast_output_frequency", 0);
    }

    const amrex::Vector<amrex::Real> m_vel{5.0_rt, 0.0_rt, 0.0_rt};
    const amrex::Real m_tau{30.0_rt};
    const amrex::Real m_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
};

TEST_F(MetMastTest, single_point)
{
    write_file("metmast_single.txt", "425 425 112.5 10 2 1 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_single.txt"));
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    EXPECT_EQ(forcing.num_stations(), 1);
    evaluate(forcing);

    // At the point
    EXPECT_NEAR(probe(8, 8, 4, 0), expected_source(0, 10, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 4, 1), expected_source(0, 2, 0, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 4, 2), expected_source(0, 1, 0, m_tau), m_tol);
    // Anisotropic Gaussian: 100 m away horizontally, 50 m vertically
    const amrex::Real ri2 = (0.2_rt * 0.2_rt) + (2.0_rt * 2.0_rt);
    EXPECT_NEAR(probe(10, 8, 6, 0), expected_source(ri2, 10, 5, m_tau), m_tol);
}

TEST_F(MetMastTest, multiple_points)
{
    write_file(
        "metmast_multiple.txt",
        "125 125 112.5 8 0 0 300\n675 675 112.5 2 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_multiple.txt"));
        pp.add("metmast_horizontal_radius", 50.0_rt);
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    EXPECT_EQ(forcing.num_stations(), 2);
    evaluate(forcing);

    // Every point is forced, not only the first one
    EXPECT_NEAR(probe(2, 2, 4, 0), expected_source(0, 8, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(13, 13, 4, 0), expected_source(0, 2, 5, m_tau), m_tol);
}

TEST_F(MetMastTest, overlapping_points)
{
    write_file(
        "metmast_overlap.txt",
        "425 425 112.5 8 0 0 300\n425 425 112.5 6 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_overlap.txt"));
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // The weights are capped at one and the targets averaged
    EXPECT_NEAR(probe(8, 8, 4, 0), expected_source(0, 7, 5, m_tau), m_tol);
}

TEST_F(MetMastTest, lidar_profile)
{
    write_file(
        "lidar_profile.txt",
        "425 425\n"
        "50 4 1 0 0 0 0\n"
        "100 6 1 0 0 0 0\n"
        "200 10 1 0 0 0 0\n");
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<std::string> files{"lidar_profile.txt"};
        pp.addarr("metmast_profile_files", files);
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    EXPECT_EQ(forcing.num_stations(), 1);
    evaluate(forcing);

    // Inside the profile: linear interpolation, full weight
    EXPECT_NEAR(probe(8, 8, 3, 0), expected_source(0, 5.5, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 5, 0), expected_source(0, 7.5, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 5, 1), expected_source(0, 1, 0, m_tau), m_tol);
    // Below the lowest gate: bottom value, tapered over 37.5 m
    EXPECT_NEAR(
        probe(8, 8, 0, 0), expected_source(1.5_rt * 1.5_rt, 4, 5, m_tau),
        m_tol);
    // Above the highest gate: top value, tapered over 62.5 m
    EXPECT_NEAR(
        probe(8, 8, 10, 0), expected_source(2.5_rt * 2.5_rt, 10, 5, m_tau),
        m_tol);
}

TEST_F(MetMastTest, sigma_band)
{
    write_file("lidar_sigma.txt", "425 425\n87.5 7 0 0 1 0 0\n");
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<std::string> files{"lidar_sigma.txt"};
        pp.addarr("metmast_profile_files", files);
    }
    setup_sim();

    // Default band of one standard deviation: forced to the band edge 6
    {
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_NEAR(probe(8, 8, 3, 0), expected_source(0, 6, 5, m_tau), m_tol);
    }
    // Velocity inside a wider band: no forcing
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_sigma_factor", 3.0_rt);
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_NEAR(probe(8, 8, 3, 0), 0.0_rt, m_tol);
    }
    // No band: plain relaxation to the mean
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_sigma_factor", 0.0_rt);
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_NEAR(probe(8, 8, 3, 0), expected_source(0, 7, 5, m_tau), m_tol);
    }
}

TEST_F(MetMastTest, terrain)
{
    write_file("metmast_terrain.txt", "425 425 62.5 10 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_terrain.txt"));
    }
    setup_sim();
    // Flat terrain 100 m high: the four lowest cells are inside it
    auto& terrain_height = sim().repo().declare_field("terrain_height", 1, 1);
    auto& terrain_blank = sim().repo().declare_int_field("terrain_blank", 1, 1);
    terrain_height.setVal(100.0_rt);
    terrain_blank.setVal(0);
    terrain_blank(0).setVal(
        1, amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 15, 3)), 1);

    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // The point is 62.5 m above the terrain, at z = 162.5 m
    EXPECT_NEAR(probe(8, 8, 6, 0), expected_source(0, 10, 5, m_tau), m_tol);
    // First cell above the terrain is 12.5 m above ground
    EXPECT_NEAR(probe(8, 8, 4, 0), expected_source(4, 10, 5, m_tau), m_tol);
    // Cells inside the terrain are not forced
    for (int k = 0; k < 4; ++k) {
        EXPECT_NEAR(probe(8, 8, k, 0), 0.0_rt, m_tol);
    }
}

TEST_F(MetMastTest, rate_limit)
{
    write_file("metmast_rate.txt", "425 425 112.5 10 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_rate.txt"));
        pp.add("meso_timescale", 1.0e-2_rt);
    }
    setup_sim();
    auto& time = sim().time();
    time.new_timestep();
    time.set_current_cfl(2.0_rt, 0.0_rt, 0.0_rt);
    EXPECT_NEAR(time.delta_t(), 0.1_rt, m_tol);

    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // tau < dt: the rate is limited to 1/dt
    EXPECT_NEAR(
        probe(8, 8, 4, 0), (10.0_rt - 5.0_rt) / 0.1_rt, 1.0e2_rt * m_tol);
}

TEST_F(MetMastTest, body_force_controller)
{
    setup_mean_gate("lidar_mean.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_output_frequency", 1);
    }
    setup_sim();
    set_velocity({5.0_rt, 1.0_rt, 0.3_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());

    // First update: the filter starts from the footprint mean, no integral
    evaluate(forcing);
    EXPECT_NEAR(forcing.footprint_velocity(0, 0), 5.0_rt, m_tol);
    EXPECT_NEAR(forcing.footprint_velocity(0, 2), 0.3_rt, m_tol);
    EXPECT_NEAR(forcing.footprint_sigma(0, 0), 0.0_rt, 1.0e2_rt * m_tol);
    EXPECT_NEAR(forcing.body_force(0, 0), 2.0_rt / m_tau, m_tol);
    EXPECT_NEAR(forcing.body_force(0, 1), 1.0_rt / m_tau, m_tol);
    // No vertical forcing by default
    EXPECT_NEAR(forcing.body_force(0, 2), 0.0_rt, m_tol);
    EXPECT_NEAR(probe(8, 8, 3, 0), 2.0_rt / m_tau, m_tol);
    EXPECT_NEAR(probe(8, 8, 3, 2), 0.0_rt, m_tol);
    EXPECT_TRUE(amrex::FileExists("post_processing/metmast_body_force.txt"));

    // Second step: filtered mean, filtered second moment and the integral
    next_step();
    const amrex::Real dt = sim().time().delta_t();
    set_velocity({6.0_rt, 1.0_rt, 0.3_rt});
    evaluate(forcing);
    const amrex::Real alpha = dt / 120.0_rt;
    const amrex::Real mean = 5.0_rt + alpha;
    const amrex::Real second = 25.0_rt + (alpha * 11.0_rt);
    const amrex::Real err = 7.0_rt - mean;
    EXPECT_NEAR(forcing.footprint_velocity(0, 0), mean, m_tol);
    EXPECT_NEAR(
        forcing.footprint_sigma(0, 0), std::sqrt(second - (mean * mean)),
        1.0e2_rt * m_tol);
    const amrex::Real force = (err / m_tau) + (err * dt / (4.0_rt * m_tau));
    EXPECT_NEAR(forcing.body_force(0, 0), force, m_tol);
    EXPECT_NEAR(probe(8, 8, 3, 0), force, m_tol);
    // Off the gate the force is spread with the station weights
    EXPECT_NEAR(
        probe(10, 8, 3, 0), force * std::exp(-0.25_rt * 0.04_rt), m_tol);
}

TEST_F(MetMastTest, body_force_independent_of_local_velocity)
{
    setup_mean_gate("lidar_local.txt");
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);
    const amrex::Real before = probe(8, 8, 3, 0);

    // A local gust in the same step does not change the force: the footprint
    // is averaged once per step and the force ignores the local velocity
    const amrex::Box gust(amrex::IntVect(7, 7, 2), amrex::IntVect(9, 9, 4));
    auto& velocity = sim().repo().get_field("velocity");
    velocity(0).setVal(50.0_rt, gust, 0, AMREX_SPACEDIM);
    velocity.state(kynema_sgf::FieldState::Old)(0).setVal(
        50.0_rt, gust, 0, AMREX_SPACEDIM);
    evaluate(forcing);
    EXPECT_NEAR(probe(8, 8, 3, 0), before, m_tol);
    EXPECT_NEAR(before, 2.0_rt / m_tau, m_tol);
}

TEST_F(MetMastTest, body_force_vertical)
{
    setup_mean_gate("lidar_vertical.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_force_vertical", true);
        pp.add("metmast_integral_timescale", 0.0_rt);
    }
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.3_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);
    next_step();
    evaluate(forcing);
    // Vertical forced, and no integral term after a step
    EXPECT_NEAR(forcing.body_force(0, 2), -0.3_rt / m_tau, m_tol);
    EXPECT_NEAR(forcing.body_force(0, 0), 2.0_rt / m_tau, m_tol);
}

TEST_F(MetMastTest, body_force_max_force)
{
    setup_mean_gate("lidar_max_force.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_max_force", 0.01_rt);
        // No time filter, so the mean follows the velocity immediately
        pp.add("metmast_averaging_time", 1.0e-3_rt);
    }
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);
    next_step();
    evaluate(forcing);
    // The force is capped
    EXPECT_NEAR(forcing.body_force(0, 0), 0.01_rt, m_tol);
    EXPECT_NEAR(probe(8, 8, 3, 0), 0.01_rt, m_tol);

    // The integral did not wind up while saturated: at the target only the
    // proportional part is left, and it is zero
    next_step();
    set_velocity({7.0_rt, 2.0_rt, 0.0_rt});
    evaluate(forcing);
    EXPECT_NEAR(forcing.body_force(0, 0), 0.0_rt, m_tol);
    EXPECT_NEAR(forcing.body_force(0, 1), 0.0_rt, m_tol);
}

TEST_F(MetMastTest, body_force_averaging_radius)
{
    setup_mean_gate("lidar_avg_radius.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_averaging_radius", 1.0_rt);
    }
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    // Faster column through the lidar only
    const amrex::Box column(amrex::IntVect(8, 8, 0), amrex::IntVect(8, 8, 15));
    auto& velocity = sim().repo().get_field("velocity");
    velocity.state(kynema_sgf::FieldState::Old)(0).setVal(9.0_rt, column, 0, 1);
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // The small average sees the lidar column only; the force keeps R_h
    EXPECT_NEAR(forcing.footprint_velocity(0, 0), 9.0_rt, m_tol);
    EXPECT_NEAR(forcing.body_force(0, 0), -2.0_rt / m_tau, m_tol);
    EXPECT_NEAR(
        probe(10, 8, 3, 0), -2.0_rt / m_tau * std::exp(-0.25_rt * 0.04_rt),
        m_tol);
}

TEST_F(MetMastTest, body_force_gate_length)
{
    setup_mean_gate("lidar_gate_length.txt");
    setup_sim();
    // Faster layer in the gate's own cells (k = 3, z = 87.5 m)
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    const amrex::Box layer(amrex::IntVect(0, 0, 3), amrex::IntVect(15, 15, 3));
    sim()
        .repo()
        .get_field("velocity")
        .state(kynema_sgf::FieldState::Old)(0)
        .setVal(9.0_rt, layer, 0, 1);

    // Default: the vertical forcing weights also average the layers around
    {
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_LT(forcing.footprint_velocity(0, 0), 8.0_rt);
    }
    // A 40 m range gate on 25 m cells averages the gate's own layer only;
    // the next layers are 25 m away
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_gate_length", 40.0_rt);
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_NEAR(forcing.footprint_velocity(0, 0), 9.0_rt, m_tol);
        EXPECT_NEAR(forcing.body_force(0, 0), -2.0_rt / m_tau, m_tol);
    }
}

TEST_F(MetMastTest, body_force_monitor)
{
    setup_mean_gate("lidar_monitor.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_forcing_type", std::string("monitor"));
    }
    setup_sim();
    set_velocity({5.0_rt, 1.0_rt, 0.0_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);
    next_step();
    evaluate(forcing);

    // The footprint average is computed, but no force is applied
    EXPECT_NEAR(forcing.footprint_velocity(0, 0), 5.0_rt, m_tol);
    EXPECT_NEAR(forcing.footprint_velocity(0, 1), 1.0_rt, m_tol);
    EXPECT_NEAR(forcing.body_force(0, 0), 0.0_rt, m_tol);
    EXPECT_NEAR(utils::field_max(src_term(), 0), 0.0_rt, m_tol);
    EXPECT_NEAR(utils::field_min(src_term(), 0), 0.0_rt, m_tol);
}

TEST_F(MetMastTest, body_force_gain_schedule)
{
    setup_mean_gate("lidar_schedule.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_gain_schedule", true);
    }
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // T_p = R_h / U = 500 / 5 s, tau_I = 2 T_p (T_avg + T_p), tau = tau_I / 10
    const amrex::Real tp = 100.0_rt;
    const amrex::Real ti = 2.0_rt * tp * (120.0_rt + tp);
    EXPECT_NEAR(forcing.level_integral_timescale(0), ti, 1.0e3_rt * m_tol);
    EXPECT_NEAR(forcing.level_timescale(0), ti / 10.0_rt, 1.0e2_rt * m_tol);
    EXPECT_NEAR(forcing.body_force(0, 0), 2.0_rt * 10.0_rt / ti, m_tol);
    EXPECT_NEAR(probe(8, 8, 3, 0), 2.0_rt * 10.0_rt / ti, m_tol);
}

TEST_F(MetMastTest, body_force_gain_schedule_floor)
{
    setup_mean_gate("lidar_schedule_floor.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_gain_schedule", true);
    }
    setup_sim();
    set_velocity({0.1_rt, 0.0_rt, 0.0_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // Nearly stagnant air: the speed is floored at 0.5 m/s
    const amrex::Real tp = 500.0_rt / 0.5_rt;
    const amrex::Real ti = 2.0_rt * tp * (120.0_rt + tp);
    EXPECT_NEAR(forcing.level_integral_timescale(0), ti, 1.0e5_rt * m_tol);
    EXPECT_NEAR(
        forcing.body_force(0, 0), (7.0_rt - 0.1_rt) * 10.0_rt / ti, m_tol);
}

TEST_F(MetMastTest, body_force_gain_schedule_integral)
{
    setup_mean_gate("lidar_schedule_integral.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_gain_schedule", true);
        // No time filter, so the gains follow the velocity immediately
        pp.add("metmast_averaging_time", 1.0e-3_rt);
    }
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);
    next_step();
    const amrex::Real dt = sim().time().delta_t();
    evaluate(forcing);
    next_step();
    set_velocity({10.0_rt, 0.0_rt, 0.0_rt});
    evaluate(forcing);

    // The integral is a force: each step adds err dt / tau_I of that step,
    // so a change of gains does not rescale the past
    const amrex::Real ti1 = 2.0_rt * 100.0_rt * (1.0e-3_rt + 100.0_rt);
    const amrex::Real ti2 = 2.0_rt * 50.0_rt * (1.0e-3_rt + 50.0_rt);
    const amrex::Real integ = (2.0_rt * dt / ti1) + (-3.0_rt * dt / ti2);
    EXPECT_NEAR(
        forcing.body_force(0, 0), (-3.0_rt * 10.0_rt / ti2) + integ, m_tol);
}

TEST_F(MetMastTest, body_force_start_time)
{
    setup_mean_gate("lidar_start_time.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_start_time", 0.25_rt);
    }
    setup_sim();
    // Spin-up transient: an empty footprint
    set_velocity({0.0_rt, 0.0_rt, 0.0_rt});
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    auto& time = sim().time();
    int nsteps = 0;
    while (time.current_time() < 0.25_rt) {
        evaluate(forcing);
        // Not started: no force and no average
        EXPECT_NEAR(forcing.body_force(0, 0), 0.0_rt, m_tol);
        EXPECT_NEAR(utils::field_max(src_term(), 0), 0.0_rt, m_tol);
        EXPECT_NEAR(forcing.footprint_velocity(0, 0), 0.0_rt, m_tol);
        next_step();
        ASSERT_LT(++nsteps, 100);
    }
    EXPECT_GT(nsteps, 0);

    // First update after the start: the filter starts from the current flow,
    // not from the transient, and the integral from zero
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    evaluate(forcing);
    EXPECT_NEAR(forcing.footprint_velocity(0, 0), 5.0_rt, m_tol);
    EXPECT_NEAR(forcing.body_force(0, 0), 2.0_rt / m_tau, m_tol);
}

TEST_F(MetMastTest, body_force_terrain)
{
    write_file("lidar_mean_terrain.txt", "425 425\n62.5 7 0 0 1 1 1\n");
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<std::string> files{"lidar_mean_terrain.txt"};
        pp.addarr("metmast_profile_files", files);
        pp.add("metmast_forcing_type", std::string("body_force"));
        pp.add("metmast_output_frequency", 0);
    }
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    // Flat terrain 100 m high with a wrong velocity inside it
    auto& terrain_height = sim().repo().declare_field("terrain_height", 1, 1);
    auto& terrain_blank = sim().repo().declare_int_field("terrain_blank", 1, 1);
    terrain_height.setVal(100.0_rt);
    terrain_blank.setVal(0);
    const amrex::Box inside(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 15, 3));
    terrain_blank(0).setVal(1, inside, 1);
    sim()
        .repo()
        .get_field("velocity")
        .state(kynema_sgf::FieldState::Old)(0)
        .setVal(100.0_rt, inside, 0, AMREX_SPACEDIM);

    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // The footprint average skips the terrain cells
    EXPECT_NEAR(forcing.footprint_velocity(0, 0), 5.0_rt, m_tol);
    // The gate is 62.5 m above the terrain
    EXPECT_NEAR(probe(8, 8, 6, 0), 2.0_rt / m_tau, m_tol);
    for (int k = 0; k < 4; ++k) {
        EXPECT_NEAR(probe(8, 8, k, 0), 0.0_rt, m_tol);
    }
}

TEST_F(MetMastTest, body_force_restart_state)
{
    write_file(
        "metmast_state_test.txt",
        "1\n"
        "6 37 10\n"
        "1 2 0\n"
        "0 0 0\n");
    setup_mean_gate("lidar_restart.txt");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_restart_state", std::string("metmast_state_test.txt"));
    }
    setup_sim();
    set_velocity({5.0_rt, 0.0_rt, 0.0_rt});
    next_step();
    const amrex::Real dt = sim().time().delta_t();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // The filter and the integral continue from the saved state
    const amrex::Real alpha = dt / 120.0_rt;
    const amrex::Real mean = 6.0_rt - alpha;
    const amrex::Real err = 7.0_rt - mean;
    EXPECT_NEAR(forcing.footprint_velocity(0, 0), mean, m_tol);
    EXPECT_NEAR(
        forcing.body_force(0, 0),
        (err / m_tau) + ((10.0_rt + (err * dt)) / (4.0_rt * m_tau)), m_tol);
}

} // namespace kynema_sgf_tests
