#include "ks_test_utils/AmrexTest.H"
#include "abl_test_utils.H"
#include "src/wind_energy/ABLWallFunction.H"
#include "src/utilities/FieldPlaneAveraging.H"
#include "src/wind_energy/MOData.H"
#include "src/wind_energy/SurfaceStability.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace ss = kynema_sgf::surface_stability;

namespace {

// Reference state shared by the pinned values. The pinned values come from a
// bracketed root of the Monin-Obukhov equations in
// tools/surface_stability_sweep.py, which is written independently of
// surface_stability::solve.
constexpr amrex::Real wspd = 8.0_rt;
constexpr amrex::Real theta = 300.0_rt;
constexpr amrex::Real zref = 15.0_rt;
constexpr amrex::Real z0 = 0.1_rt;

constexpr amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
// solve() stops once utau changes by less than 1e-5 m/s; the larger bounds
// apply to single precision
constexpr amrex::Real rel_tol =
    amrex::max<amrex::Real>(1.0e-4_rt, 1.0e3_rt * eps);
//! Absolute tolerance on temperatures of about 300 K
constexpr amrex::Real temp_tol =
    amrex::max<amrex::Real>(1.0e-4_rt, 2.0e5_rt * eps);

ss::StabilityParams make_params(const ss::StabilityType type)
{
    ss::StabilityParams p;
    p.type = type;
    p.non_neutral = true;
    p.kappa = 0.41_rt;
    p.gravity = 9.81_rt;
    return p;
}

void expect_rel(const amrex::Real value, const amrex::Real expected)
{
    EXPECT_NEAR(value, expected, rel_tol * std::abs(expected));
}

class SurfaceStabilityTest : public AmrexTest
{};

} // namespace

TEST(SurfaceStability, psi_pinned)
{
    // Dyer (1974) functions with beta = 16 and gamma = 5
    constexpr amrex::Real tol = 100.0_rt * eps;
    EXPECT_NEAR(ss::psi_m(-1.0_rt, 16.0_rt, 5.0_rt), 1.11623224976833_rt, tol);
    EXPECT_NEAR(ss::psi_h(-1.0_rt, 16.0_rt, 5.0_rt), 1.88122728421442_rt, tol);
    EXPECT_NEAR(ss::psi_m(-0.1_rt, 16.0_rt, 5.0_rt), 0.283613711212781_rt, tol);
    EXPECT_NEAR(ss::psi_h(-0.1_rt, 16.0_rt, 5.0_rt), 0.534283781948425_rt, tol);
    EXPECT_NEAR(ss::psi_m(0.0_rt, 16.0_rt, 5.0_rt), 0.0_rt, tol);
    EXPECT_NEAR(ss::psi_h(0.0_rt, 16.0_rt, 5.0_rt), 0.0_rt, tol);
    EXPECT_NEAR(ss::psi_m(0.5_rt, 16.0_rt, 5.0_rt), -2.5_rt, tol);
    EXPECT_NEAR(ss::psi_h(0.5_rt, 16.0_rt, 5.0_rt), -2.5_rt, tol);

    // MOData uses the same functions
    for (const amrex::Real zeta : {-2.0_rt, -0.3_rt, 0.0_rt, 0.7_rt}) {
        EXPECT_EQ(
            kynema_sgf::MOData::calc_psi_m(zeta, 16.0_rt, 5.0_rt),
            ss::psi_m(zeta, 16.0_rt, 5.0_rt));
        EXPECT_EQ(
            kynema_sgf::MOData::calc_psi_h(zeta, 16.0_rt, 5.0_rt),
            ss::psi_h(zeta, 16.0_rt, 5.0_rt));
    }
}

TEST(SurfaceStability, heat_flux_unstable)
{
    auto p = make_params(ss::StabilityType::HEAT_FLUX);
    p.surf_temp_flux = 0.05_rt;
    const auto s = ss::solve(p, wspd, theta, zref, z0);
    expect_rel(s.utau, 0.670103665964_rt);
    expect_rel(s.obukhov_len, -448.873914456_rt);
    EXPECT_NEAR(s.surf_temp, 300.870821415_rt, temp_tol);
    EXPECT_DOUBLE_EQ(s.surf_temp_flux, 0.05_rt);
}

TEST(SurfaceStability, heat_flux_stable)
{
    auto p = make_params(ss::StabilityType::HEAT_FLUX);
    p.surf_temp_flux = -0.01_rt;
    const auto s = ss::solve(p, wspd, theta, zref, z0);
    expect_rel(s.utau, 0.649855723245_rt);
    expect_rel(s.obukhov_len, 2047.00656225_rt);
    EXPECT_NEAR(s.surf_temp, 299.810566802_rt, temp_tol);
}

TEST(SurfaceStability, heat_flux_low_wind_unstable)
{
    // The root zeta = -13.3 lies below the lower bound, where
    // ln(zref/z0) - psi_h = ln(zref/z0)/2 (zeta = -2.187). zeta is held at the
    // bound and the heat flux is reduced to what the bound carries. Pressel
    // review case F2 (U = 0.5 m/s) used to converge to the u* floor with a
    // surface temperature of the wrong sign.
    auto p = make_params(ss::StabilityType::HEAT_FLUX);
    p.surf_temp_flux = 0.4_rt;
    const auto s = ss::solve(p, 1.0_rt, theta, zref, z0);
    expect_rel(zref / s.obukhov_len, -2.18695330058_rt);
    expect_rel(s.utau, 0.118378563375_rt);
    expect_rel(s.surf_temp_flux, 0.0180399710526_rt);
    EXPECT_NEAR(s.surf_temp, 300.931197308_rt, temp_tol);

    p.surf_temp_flux = 0.1_rt;
    const auto f2 = ss::solve(p, 0.5_rt, theta, 18.75_rt, z0);
    expect_rel(f2.utau, 0.0568337509211_rt);
    EXPECT_GT(f2.surf_temp, theta);
    EXPECT_NEAR(f2.surf_temp, 300.204881127_rt, temp_tol);
}

TEST(SurfaceStability, heat_flux_stable_two_roots)
{
    // Both roots of the stable heat-flux residual lie in [0, zeta_max]; the
    // search must return the first (physical) one, zeta = 0.196, and not
    // treat the surface as decoupled
    auto p = make_params(ss::StabilityType::HEAT_FLUX);
    p.surf_temp_flux = -0.02_rt;
    const auto s = ss::solve(p, 4.0_rt, theta, zref, z0);
    expect_rel(s.utau, 0.273743621993_rt);
    expect_rel(s.obukhov_len, 76.5015831549_rt);
    EXPECT_DOUBLE_EQ(s.surf_temp_flux, -0.02_rt);
    EXPECT_NEAR(s.surf_temp, 298.932416216_rt, temp_tol);
}

TEST(SurfaceStability, heat_flux_below_critical_wind)
{
    // Pressel review case F1: a stable flux of -0.02 K m/s has no
    // Monin-Obukhov solution at 3 m/s and 18.75 m. zeta is held at
    // zeta_max = 2 and the flux is reduced, so the surface temperature stays
    // within a few kelvin instead of -1e13 K.
    auto p = make_params(ss::StabilityType::HEAT_FLUX);
    p.surf_temp_flux = -0.02_rt;
    const auto s = ss::solve(p, 3.0_rt, theta, 18.75_rt, z0);
    expect_rel(s.obukhov_len, 18.75_rt / ss::zeta_max);
    expect_rel(s.utau, 0.0807416211356_rt);
    expect_rel(s.surf_temp_flux, -0.0041878344877_rt);
    EXPECT_NEAR(s.surf_temp, 298.072848604_rt, temp_tol);
}

TEST(SurfaceStability, surface_temperature_decoupled)
{
    // Beyond the critical bulk Richardson number zeta is held at zeta_max and
    // the flux follows from the surface temperature
    auto p = make_params(ss::StabilityType::SURFACE_TEMPERATURE);
    p.surf_temp = 290.0_rt;
    const auto s = ss::solve(p, 1.0_rt, theta, 18.75_rt, z0);
    expect_rel(s.obukhov_len, 18.75_rt / ss::zeta_max);
    expect_rel(s.utau, 0.0269138737119_rt);
    expect_rel(s.surf_temp_flux, -0.00724356598179_rt);
}

TEST(SurfaceStability, surface_temperature)
{
    for (const auto type :
         {ss::StabilityType::SURFACE_TEMPERATURE,
          ss::StabilityType::HEATING_RATE}) {
        auto p = make_params(type);
        p.surf_temp = 302.0_rt;
        const auto su = ss::solve(p, wspd, theta, zref, z0);
        expect_rel(su.utau, 0.686271729429_rt);
        expect_rel(su.obukhov_len, -195.814871651_rt);
        expect_rel(su.surf_temp_flux, 0.123115021299_rt);
        EXPECT_DOUBLE_EQ(su.surf_temp, 302.0_rt);

        p.surf_temp = 298.0_rt;
        const auto st = ss::solve(p, wspd, theta, zref, z0);
        expect_rel(st.utau, 0.604438074663_rt);
        expect_rel(st.obukhov_len, 180.335071131_rt);
        expect_rel(st.surf_temp_flux, -0.0913363465255_rt);
    }
}

TEST(SurfaceStability, obukhov_length_round_trip)
{
    // The Obukhov length of the heat-flux case gives back its heat flux and
    // surface temperature in closed form
    auto p = make_params(ss::StabilityType::MONIN_OBUKHOV_LENGTH);
    p.obukhov_len = -448.873914456_rt;
    const auto s = ss::solve(p, wspd, theta, zref, z0);
    expect_rel(s.utau, 0.670103665964_rt);
    expect_rel(s.surf_temp_flux, 0.05_rt);
    EXPECT_NEAR(s.surf_temp, 300.870821415_rt, temp_tol);
}

TEST(SurfaceStability, neutral)
{
    const amrex::Real utau_neutral = 0.41_rt * wspd / std::log(zref / z0);
    {
        auto p = make_params(ss::StabilityType::HEAT_FLUX);
        p.surf_temp_flux = 0.0_rt;
        const auto s = ss::solve(p, wspd, theta, zref, z0);
        EXPECT_DOUBLE_EQ(s.utau, utau_neutral);
        EXPECT_DOUBLE_EQ(s.obukhov_len, kynema_sgf::constants::LARGE_NUM);
        EXPECT_DOUBLE_EQ(s.surf_temp, theta);
    }
    {
        auto p = make_params(ss::StabilityType::SURFACE_TEMPERATURE);
        p.surf_temp = theta;
        const auto s = ss::solve(p, wspd, theta, zref, z0);
        EXPECT_DOUBLE_EQ(s.utau, utau_neutral);
        EXPECT_DOUBLE_EQ(s.surf_temp_flux, 0.0_rt);
    }
    {
        // wall_het_model is not mol
        auto p = make_params(ss::StabilityType::HEAT_FLUX);
        p.non_neutral = false;
        p.surf_temp_flux = 0.2_rt;
        const auto s = ss::solve(p, wspd, theta, zref, z0);
        EXPECT_DOUBLE_EQ(s.utau, utau_neutral);
    }
}

TEST(SurfaceStability, profile_passes_reference)
{
    for (const amrex::Real q : {-0.01_rt, 0.05_rt, 0.3_rt}) {
        auto p = make_params(ss::StabilityType::HEAT_FLUX);
        p.surf_temp_flux = q;
        const auto s = ss::solve(p, wspd, theta, zref, z0);
        EXPECT_NEAR(ss::theta_at(p, s, zref, z0), theta, 10.0_rt * temp_tol);
    }
}

TEST(SurfaceStability, zero_wind_is_finite)
{
    for (const auto type :
         {ss::StabilityType::HEAT_FLUX,
          ss::StabilityType::SURFACE_TEMPERATURE}) {
        auto p = make_params(type);
        p.surf_temp_flux = 0.1_rt;
        p.surf_temp = 305.0_rt;
        const auto s = ss::solve(p, 0.0_rt, theta, zref, z0);
        EXPECT_TRUE(std::isfinite(s.utau));
        EXPECT_TRUE(std::isfinite(s.obukhov_len));
        EXPECT_TRUE(std::isfinite(s.surf_temp));
        EXPECT_TRUE(std::isfinite(s.surf_temp_flux));
    }
}

TEST_F(SurfaceStabilityTest, parse_default)
{
    const auto p = ss::parse("test", 290.0_rt);
    EXPECT_EQ(p.type, ss::StabilityType::MONIN_OBUKHOV_LENGTH);
    EXPECT_FALSE(p.non_neutral);
    EXPECT_FALSE(p.is_local());
    EXPECT_DOUBLE_EQ(p.obukhov_len, kynema_sgf::constants::LARGE_NUM);
}

TEST_F(SurfaceStabilityTest, parse_heat_flux)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("heat_flux"));
    pp.add("surface_temp_flux", 0.05_rt);
    pp.add("mo_gamma_h", 4.0_rt);
    pp.add("mo_beta_h", 12.0_rt);
    const auto p = ss::parse("test", 290.0_rt);
    EXPECT_EQ(p.type, ss::StabilityType::HEAT_FLUX);
    EXPECT_TRUE(p.is_local());
    EXPECT_DOUBLE_EQ(p.surf_temp_flux, 0.05_rt);
    // The heat constants are read from their own keys
    EXPECT_DOUBLE_EQ(p.gamma_h, 4.0_rt);
    EXPECT_DOUBLE_EQ(p.beta_h, 12.0_rt);
    EXPECT_DOUBLE_EQ(p.gamma_m, 5.0_rt);
    EXPECT_DOUBLE_EQ(p.beta_m, 16.0_rt);
}

TEST_F(SurfaceStabilityTest, parse_heating_rate)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("heating_rate"));
    pp.add("surface_temp_rate", -0.5_rt);
    pp.add("surface_temp_rate_tstart", 3600.0_rt);
    auto p = ss::parse("test", 290.0_rt);
    EXPECT_EQ(p.type, ss::StabilityType::HEATING_RATE);
    // surface_temp_init defaults to the reference temperature
    EXPECT_DOUBLE_EQ(p.surf_temp, 290.0_rt);
    ss::update(p, 1800.0_rt);
    EXPECT_DOUBLE_EQ(p.surf_temp, 290.0_rt);
    ss::update(p, 3600.0_rt + 7200.0_rt);
    EXPECT_DOUBLE_EQ(p.surf_temp, 289.0_rt);
}

TEST_F(SurfaceStabilityTest, parse_surface_temperature)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("surface_temperature"));
    pp.add("surface_temp_init", 302.0_rt);
    auto p = ss::parse("test", 290.0_rt);
    EXPECT_EQ(p.type, ss::StabilityType::SURFACE_TEMPERATURE);
    ss::update(p, 1.0e5_rt);
    EXPECT_DOUBLE_EQ(p.surf_temp, 302.0_rt);
}

TEST_F(SurfaceStabilityTest, parse_rejects_unknown_type)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("heatflux"));
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_needs_mol)
{
    amrex::ParmParse pp("ABL");
    pp.add("surface_stability_type", std::string("heat_flux"));
    pp.add("surface_temp_flux", 0.05_rt);
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_rejects_zero_length)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("monin_obukhov_length", 0.0_rt);
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_heat_flux_range)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("heat_flux"));
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
    pp.add("surface_temp_flux", 11.0_rt);
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_surface_temperature_checks)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("surface_temperature"));
    // surface_temp_init is required
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
    pp.add("surface_temp_init", -1.0_rt);
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_surface_temperature_rejects_rate)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("surface_temperature"));
    pp.add("surface_temp_init", 302.0_rt);
    pp.add("surface_temp_rate", 0.5_rt);
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_rejects_conflicting_keys)
{
    // The flat-wall ABLWallFunction would use the heat flux while the terrain
    // used the heating rate
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("heating_rate"));
    pp.add("surface_temp_rate", -0.25_rt);
    pp.add("surface_temp_flux", 0.05_rt);
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_heat_flux_rejects_timetable)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("heat_flux"));
    pp.add("surface_temp_flux", 0.05_rt);
    pp.add("surface_temp_flux_timetable", std::string("flux.txt"));
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, parse_zero_length_without_mol)
{
    // L is not used without wall_het_model = mol
    amrex::ParmParse pp("ABL");
    pp.add("monin_obukhov_length", 0.0_rt);
    EXPECT_NO_THROW(ss::parse("test", 300.0_rt));
}

TEST(SurfaceStability, roughness_above_reference)
{
    // z0 >= zref is limited to zref/e and gives finite fluxes of the right sign
    auto p = make_params(ss::StabilityType::HEAT_FLUX);
    p.surf_temp_flux = 0.05_rt;
    const auto s = ss::solve(p, wspd, theta, 1.5_rt, 2.0_rt);
    EXPECT_TRUE(std::isfinite(s.utau));
    EXPECT_GT(s.utau, 0.0_rt);
    EXPECT_LT(s.obukhov_len, 0.0_rt);
    EXPECT_GT(s.surf_temp, theta);
}

TEST_F(SurfaceStabilityTest, parse_heating_rate_tstart)
{
    amrex::ParmParse pp("ABL");
    pp.add("wall_het_model", std::string("mol"));
    pp.add("surface_stability_type", std::string("heating_rate"));
    pp.add("surface_temp_rate", 0.5_rt);
    pp.add("surface_temp_rate_tstart", -1.0_rt);
    EXPECT_THROW(ss::parse("test", 300.0_rt), amrex::RuntimeError);
}

TEST_F(SurfaceStabilityTest, require_obukhov_length)
{
    auto p = make_params(ss::StabilityType::MONIN_OBUKHOV_LENGTH);
    EXPECT_NO_THROW(ss::require_obukhov_length(p, "test"));
    for (const auto type :
         {ss::StabilityType::HEAT_FLUX, ss::StabilityType::HEATING_RATE,
          ss::StabilityType::SURFACE_TEMPERATURE}) {
        p.type = type;
        EXPECT_THROW(
            ss::require_obukhov_length(p, "test"), amrex::RuntimeError);
    }
}

/** Compare the per-cell solve with the flat-wall (non-immersed-boundary)
 *  ABLWallFunction, which solves the same surface layer from plane averages
 *  with MOData::update_fluxes. Both read the same ABL inputs.
 */
class SurfaceStabilityWallFunctionTest : public ABLMeshTest
{
protected:
    void populate_parameters() override
    {
        ABLMeshTest::populate_parameters();
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 2.5_rt);
        amrex::ParmParse pp_abl("ABL");
        pp_abl.add("wall_het_model", std::string("mol"));
    }

    //! Uniform 8 m/s wind at 300 K
    void setup_fields()
    {
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        pde_mgr.register_transport_pde("Temperature");
        auto& velocity = sim().repo().get_field("velocity");
        velocity.setVal(0.0_rt);
        velocity.setVal(8.0_rt, 0, 1);
        sim().repo().get_field("temperature").setVal(300.0_rt);
    }

    void advance_one_step()
    {
        auto& time = sim().time();
        time.new_timestep();
        time.set_current_cfl(0.0_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        time.new_timestep();
    }

    //! Surface layer of the wall function at the current time
    kynema_sgf::MOData wall_function_state()
    {
        kynema_sgf::VelPlaneAveraging vpa(sim(), 2, -1);
        kynema_sgf::FieldPlaneAveraging tpa(
            sim().repo().get_field("temperature"), sim().time(), 2, -1);
        vpa();
        tpa();
        kynema_sgf::ABLWallFunction wall_func(sim());
        wall_func.init_log_law_height(0);
        wall_func.update_umean(vpa, tpa);
        return wall_func.mo();
    }

    //! Heat-flux input read by both the wall function and parse()
    void check_heat_flux(const amrex::Real q)
    {
        populate_parameters();
        {
            amrex::ParmParse pp("ABL");
            pp.add("surface_stability_type", std::string("heat_flux"));
            pp.add("surface_temp_flux", q);
        }
        setup_fields();
        const auto mo = wall_function_state();
        EXPECT_EQ(mo.alg_type, kynema_sgf::MOData::ThetaCalcType::HEAT_FLUX);
        const auto p = ss::parse("test", 300.0_rt);
        EXPECT_EQ(p.type, ss::StabilityType::HEAT_FLUX);
        compare(p, mo);
    }

    //! Solve the same surface layer with surface_stability::solve
    static void
    compare(const ss::StabilityParams& p, const kynema_sgf::MOData& mo)
    {
        // MOData::update_fluxes and solve() both stop once utau changes by
        // less than 1e-5 m/s; MOData also lags psi_h by one iteration
        constexpr amrex::Real tol = 1.0e-3_rt;
        EXPECT_DOUBLE_EQ(mo.z0t, mo.z0);
        EXPECT_DOUBLE_EQ(mo.alpha_h, 1.0_rt);
        const auto s =
            ss::solve(p, mo.vmag_mean, mo.theta_mean, mo.zref, mo.z0);
        EXPECT_NEAR(s.utau, mo.utau, tol * mo.utau);
        EXPECT_NEAR(
            s.obukhov_len, mo.obukhov_len, tol * std::abs(mo.obukhov_len));
        EXPECT_NEAR(
            s.surf_temp_flux, mo.surf_temp_flux,
            tol * std::abs(mo.surf_temp_flux));
        EXPECT_NEAR(s.surf_temp, mo.surf_temp, tol);
    }
};

TEST_F(SurfaceStabilityWallFunctionTest, heat_flux_unstable)
{
    check_heat_flux(0.05_rt);
}

TEST_F(SurfaceStabilityWallFunctionTest, heat_flux_stable)
{
    check_heat_flux(-0.01_rt);
}

TEST_F(SurfaceStabilityWallFunctionTest, heat_flux_strongly_unstable)
{
    check_heat_flux(0.2_rt);
}

TEST_F(SurfaceStabilityWallFunctionTest, heating_rate)
{
    populate_parameters();
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_stability_type", std::string("heating_rate"));
        pp.add("surface_temp_rate", -3600.0_rt);
        pp.add("surface_temp_init", 300.5_rt);
        pp.add("surface_temp_rate_tstart", 0.5_rt);
    }
    setup_fields();
    // At t = 2.5 s the surface has cooled by 2 K to 298.5 K
    advance_one_step();
    const auto mo = wall_function_state();
    EXPECT_EQ(
        mo.alg_type, kynema_sgf::MOData::ThetaCalcType::SURFACE_TEMPERATURE);
    EXPECT_NEAR(mo.surf_temp, 298.5_rt, 1.0e-10_rt);
    auto p = ss::parse("test", 300.0_rt);
    ss::update(p, sim().time().current_time());
    EXPECT_NEAR(p.surf_temp, 298.5_rt, 1.0e-10_rt);
    compare(p, mo);
}

TEST_F(SurfaceStabilityWallFunctionTest, surface_temperature)
{
    // The wall function has no constant surface temperature input; a zero
    // heating rate from surface_temp_init gives the same surface
    populate_parameters();
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_rate", 0.0_rt);
        pp.add("surface_temp_init", 302.0_rt);
    }
    setup_fields();
    const auto mo = wall_function_state();
    {
        amrex::ParmParse pp("ABL");
        pp.remove("surface_temp_rate");
        pp.add("surface_stability_type", std::string("surface_temperature"));
    }
    const auto p = ss::parse("test", 300.0_rt);
    EXPECT_EQ(p.type, ss::StabilityType::SURFACE_TEMPERATURE);
    compare(p, mo);
}

TEST_F(SurfaceStabilityWallFunctionTest, obukhov_length)
{
    // Feed the Obukhov length found by the wall function back as the input
    populate_parameters();
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_flux", 0.05_rt);
    }
    setup_fields();
    const auto mo = wall_function_state();
    {
        amrex::ParmParse pp("ABL");
        pp.add("monin_obukhov_length", mo.obukhov_len);
    }
    const auto p = ss::parse("test", 300.0_rt);
    EXPECT_EQ(p.type, ss::StabilityType::MONIN_OBUKHOV_LENGTH);
    compare(p, mo);
}

} // namespace kynema_sgf_tests
