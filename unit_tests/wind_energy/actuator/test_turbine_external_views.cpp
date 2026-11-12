#include <array>
#include <vector>
#include "ks_test_utils/AmrexTest.H"
#include "test_ext_turbine_utils.H"
#include "src/wind_energy/actuator/turbine/external/turbine_external_ops.H"
#include "AMReX_Exception.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

// Stand-ins for the turbine and solver data of an external solver. The
// turbine holds the buffers exchanged with the solver: the fluid velocity
// at the velocity points and the velocity of the force points.
struct MockExtTurbine
{
    float* fluid_velocity(const int dir) { return fluid_vel[dir].data(); }
    float* solid_velocity(const int dir) { return solid_vel[dir].data(); }
    int length_fluid_velocity(const int dir) const
    {
        return static_cast<int>(fluid_vel[dir].size());
    }

    std::array<std::vector<float>, 3> fluid_vel;
    std::array<std::vector<float>, 3> solid_vel;
};

struct MockExtSolverData
{};

using MockViewData = ExtTurbineViewData<MockExtTurbine, MockExtSolverData>;

} // namespace

TEST(TurbineExternalViews, hub_first_velocity_offset_one)
{
    MockViewData data;
    init_hub_first_points(data, 3, 5, 5);

    // The hub is first in the velocity points: skip it, as in the force
    // points
    data.meta().velocity_points_offset = 1;
    kynema_sgf::actuator::external::make_component_views(data);

    check_blade_views_aligned(data);
}

TEST(TurbineExternalViews, hub_first_velocity_offset_zero)
{
    namespace vs = kynema_sgf::vs;
    const int num_blades = 3;
    const int npb = 5;
    MockViewData data;
    init_hub_first_points(data, num_blades, npb, 5);

    // Without the offset every blade view starts one velocity point early
    data.meta().velocity_points_offset = 0;
    kynema_sgf::actuator::external::make_component_views(data);

    const auto& grid = data.grid();
    const auto& hub = data.meta().hub;
    const auto& blades = data.meta().blades;
    ASSERT_EQ(static_cast<int>(blades.size()), num_blades);
    for (int ib = 0; ib < num_blades; ++ib) {
        const auto& cv = blades[ib];
        const int start = ib * npb;
        EXPECT_EQ(std::distance(grid.vel.data(), cv.vel.data()), start);
        EXPECT_EQ(
            std::distance(grid.vel_pos.data(), cv.vel_pos.data()), start);
        // Entry k is the velocity point of force point k - 1
        for (int k = 1; k < npb; ++k) {
            EXPECT_EQ(vs::mag(cv.vel_pos[k] - cv.pos[k - 1]), 0.0_rt)
                << "blade " << ib << ", node " << k;
        }
    }
    // Entry 0 is the hub for the first blade and the tip of the previous
    // blade for the other two
    EXPECT_EQ(vs::mag(blades[0].vel_pos[0] - hub.pos[0]), 0.0_rt);
    EXPECT_EQ(vs::mag(blades[1].vel_pos[0] - blades[0].pos[npb - 1]), 0.0_rt);
    EXPECT_EQ(vs::mag(blades[2].vel_pos[0] - blades[1].pos[npb - 1]), 0.0_rt);
    // The velocity point at the tip of the last blade is in no blade view
    EXPECT_EQ(
        std::distance(grid.vel.data(), blades[2].vel.end()), num_blades * npb);
}

TEST(TurbineExternalViews, fllc_relative_velocity_with_different_point_counts)
{
    namespace vs = kynema_sgf::vs;
    const int num_blades = 3;
    const int npb = 4; // force points on a blade
    const int nvb = 6; // velocity points on a blade
    MockViewData data;
    init_hub_first_points(data, num_blades, npb, nvb, 2);
    auto& grid = data.grid();
    auto& meta = data.meta();
    meta.velocity_points_offset = 1;
    kynema_sgf::actuator::external::make_component_views(data);
    data.info().is_root_proc = true;

    // Fluid velocity last sent to the external solver: linear in space, so
    // that its interpolation along a blade is exact. Velocity of the force
    // points: different for every point. Velocity sampled now: not linear.
    const auto fluid_u = [](const vs::Vector& x) {
        return 1.0_rt + 2.0_rt * x.y() + 3.0_rt * x.z();
    };
    const int npts = static_cast<int>(grid.pos.size());
    const int npts_vel = static_cast<int>(grid.vel_pos.size());
    auto& ext = meta.ext_data;
    for (int dir = 0; dir < 3; ++dir) {
        ext.fluid_vel[dir].assign(npts_vel, 0.0F);
        ext.solid_vel[dir].assign(npts, 0.0F);
    }
    for (int i = 0; i < npts_vel; ++i) {
        ext.fluid_vel[0][i] = static_cast<float>(fluid_u(grid.vel_pos[i]));
        grid.vel[i] = vs::Vector{
            static_cast<amrex::Real>(7 + (i * i) % 5), 0.0_rt, 0.0_rt};
    }
    for (int i = 0; i < npts; ++i) {
        ext.solid_vel[0][i] = 0.25F * static_cast<float>(i);
    }

    // Correction initialized and no force on the blades
    meta.fllc.resize(num_blades);
    for (int ib = 0; ib < num_blades; ++ib) {
        meta.fllc[ib].epsilon = 2.0_rt;
        kynema_sgf::actuator::fllc_init(
            meta.fllc[ib], meta.blades[ib], 1.0_rt);
    }

    kynema_sgf::actuator::ops::update_vel_op(data);

    // Relative velocity at force point k of every blade: the fluid velocity
    // at that point minus the velocity of that point
    const amrex::Real tol = 1.0e-5_rt;
    for (int ib = 0; ib < num_blades; ++ib) {
        for (int k = 0; k < npb; ++k) {
            const int ip = 1 + ib * npb + k;
            const auto gold =
                fluid_u(grid.pos[ip]) - 0.25_rt * static_cast<amrex::Real>(ip);
            EXPECT_NEAR(meta.blades[ib].vel_rel[k].x(), gold, tol)
                << "blade " << ib << ", force node " << k;
            EXPECT_EQ(meta.blades[ib].vel_rel[k].y(), 0.0_rt);
            EXPECT_EQ(meta.blades[ib].vel_rel[k].z(), 0.0_rt);
        }
    }
    // Without force there is no correction: the velocity sent to the
    // external solver is the sampled one, at every velocity point
    for (int i = 0; i < npts_vel; ++i) {
        EXPECT_EQ(ext.fluid_vel[0][i], static_cast<float>(7 + (i * i) % 5))
            << "velocity point " << i;
        EXPECT_EQ(ext.fluid_vel[1][i], 0.0F);
        EXPECT_EQ(ext.fluid_vel[2][i], 0.0F);
    }
}

TEST(TurbineExternalViews, epsilon_must_be_positive_on_blades)
{
    MockViewData data;
    init_hub_first_points(data, 3, 5, 0);
    data.meta().velocity_points_offset = 1;
    kynema_sgf::actuator::external::make_component_views(data);

    // None of the epsilon inputs given: the Gaussian of every blade point
    // has no width, and its force would not reach the mesh
    EXPECT_THROW(
        kynema_sgf::actuator::external::init_epsilon(data),
        amrex::RuntimeError);
}

TEST(TurbineExternalViews, epsilon_must_be_positive_on_tower)
{
    namespace vs = kynema_sgf::vs;
    MockViewData data;
    init_hub_first_points(data, 3, 5, 4);
    data.meta().velocity_points_offset = 1;
    kynema_sgf::actuator::external::make_component_views(data);

    // Tower points, and an epsilon for the blades only
    data.meta().eps_chord = vs::Vector{0.25_rt, 0.25_rt, 0.25_rt};
    EXPECT_THROW(
        kynema_sgf::actuator::external::init_epsilon(data),
        amrex::RuntimeError);
}

TEST(TurbineExternalViews, epsilon_from_chord_and_tower_inputs)
{
    namespace vs = kynema_sgf::vs;
    const int num_blades = 3;
    const int npb = 5;
    const int npt = 4;
    MockViewData data;
    init_hub_first_points(data, num_blades, npb, npt);
    auto& meta = data.meta();
    meta.velocity_points_offset = 1;
    kynema_sgf::actuator::external::make_component_views(data);

    meta.eps_chord = vs::Vector{0.25_rt, 0.25_rt, 0.25_rt};
    meta.eps_tower = vs::Vector{2.0_rt, 2.0_rt, 2.0_rt};
    kynema_sgf::actuator::external::init_epsilon(data);

    // Unit chord on the blades
    for (int ib = 0; ib < num_blades; ++ib) {
        for (int k = 0; k < npb; ++k) {
            EXPECT_EQ(meta.blades[ib].epsilon[k].x(), 0.25_rt);
        }
    }
    for (int j = 0; j < npt; ++j) {
        EXPECT_EQ(meta.tower.epsilon[j].x(), 2.0_rt);
    }
}

TEST(TurbineExternalViews, bound_box_holds_the_gaussians)
{
    namespace vs = kynema_sgf::vs;
    MockViewData data;
    init_hub_first_points(data, 3, 5, 4);
    auto& meta = data.meta();
    meta.velocity_points_offset = 1;
    kynema_sgf::actuator::external::make_component_views(data);

    // Box made from the inputs of a turbine at the origin with a rotor
    // radius of 5: 1.25 radii around it. The points are up to 5 away from
    // the hub
    data.info().bound_box = amrex::RealBox(
        -6.25_rt, -6.25_rt, -6.25_rt, 6.25_rt, 6.25_rt, 6.25_rt);
    meta.eps_inp = vs::Vector{2.0_rt, 2.0_rt, 2.0_rt};
    kynema_sgf::actuator::external::init_epsilon(data);

    // The Gaussian of a point is cut off at 4 epsilon. Patches outside of
    // the box receive no force, so the box has to reach that far
    const amrex::Real reach = 4.0_rt * 2.0_rt;
    const auto& bbox = data.info().bound_box;
    for (const auto& pos : data.grid().pos) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_LE(bbox.lo(n), pos[n] - reach);
            EXPECT_GE(bbox.hi(n), pos[n] + reach);
        }
    }
}

} // namespace kynema_sgf_tests
