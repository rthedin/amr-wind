#include "ks_test_utils/AmrexTest.H"
#include "test_ext_turbine_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

// Stand-ins for the turbine and solver data of an external solver
struct MockExtTurbine
{};

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

} // namespace kynema_sgf_tests
