#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/turbulence/TurbulenceModel.H"
#include "ks_test_utils/test_utils.H"
#include "src/utilities/math_ops.H"
#include "src/equation_systems/tke/source_terms/KransAxell.H"
#include "src/utilities/constants.H"
#include <array>
#include <fstream>
#include <type_traits>

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void init_strain_field(kynema_sgf::Field& fld, amrex::Real srate)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    amrex::Real offset =
        (fld.field_location() == kynema_sgf::FieldLoc::CELL) ? 0.5_rt : 0.0_rt;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + offset) * dx[0]);
                const amrex::Real y = problo[1] + ((j + offset) * dx[1]);
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);

                farrs[nbx](i, j, k, 0) = x / std::sqrt(6.0_rt) * srate;
                farrs[nbx](i, j, k, 1) = y / std::sqrt(6.0_rt) * srate;
                farrs[nbx](i, j, k, 2) = z / std::sqrt(6.0_rt) * srate;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_temperature_field(kynema_sgf::Field& fld, amrex::Real tgrad)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset =
        (fld.field_location() == kynema_sgf::FieldLoc::CELL) ? 0.5_rt : 0.0_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);

                farrs[nbx](i, j, k, 0) = z * tgrad;
            });
    }
    amrex::Gpu::streamSynchronize();
}

// value below ztop, zero above
void set_below(
    kynema_sgf::Field& fld, const amrex::Real value, const amrex::Real ztop)
{
    const auto& mesh = fld.repo().mesh();
    for (int lev = 0; lev < fld.repo().num_active_levels(); ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                farrs[nbx](i, j, k) = (z < ztop) ? value : 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}

// Neutral KLAxell inputs of test_1eqKrans_setup_calc, plus the 1-D RANS
// profile needed by the KransAxell source
void add_klaxell_inputs(
    const amrex::Real Tref, const amrex::Real gravz, const amrex::Real rho0)
{
    const std::string rans_file = "rans_1d_canopy.info";
    {
        std::ofstream os(rans_file);
        os << "0 8 0 0 0.1\n";
        os << "500 8 0 0 0.1\n";
        os << "4000 8 0 0 0.1\n";
    }
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "KLAxell");
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -gravz};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_rate", 0.0_rt);
        pp.add("initial_wind_profile", true);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 4000.0_rt};
        pp.addarr("temperature_heights", t_hts);
        pp.addarr("wind_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{Tref, Tref, Tref};
        pp.addarr("temperature_values", t_vals);
        amrex::Vector<amrex::Real> u_vals{8.0_rt, 8.0_rt, 8.0_rt};
        pp.addarr("u_values", u_vals);
        amrex::Vector<amrex::Real> v_vals{0.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("v_values", v_vals);
        amrex::Vector<amrex::Real> tke_vals{0.1_rt, 0.1_rt, 0.1_rt};
        pp.addarr("tke_values", tke_vals);
        pp.add("surface_temp_flux", 0.0_rt);
        pp.add("rans_1dprofile_file", rans_file);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }
}

} // namespace

class TurbRANSTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{10, 10, 64}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{
                {1024.0_rt, 1024.0_rt, 1024.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

TEST_F(TurbRANSTest, test_1eqKrans_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real Tref = 265.0_rt;
    const amrex::Real gravz = 10.0_rt;
    const amrex::Real rho0 = 1.2_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "KLAxell");
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -gravz};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_rate", 0.0_rt);
        pp.add("initial_wind_profile", true);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 4000.0_rt};
        pp.addarr("temperature_heights", t_hts);
        pp.addarr("wind_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{265.0_rt, 265.0_rt, 265.0_rt};
        pp.addarr("temperature_values", t_vals);
        amrex::Vector<amrex::Real> u_vals{8.0_rt, 8.0_rt, 8.0_rt};
        pp.addarr("u_values", u_vals);
        amrex::Vector<amrex::Real> v_vals{0.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("v_values", v_vals);
        amrex::Vector<amrex::Real> tke_vals{0.1_rt, 0.1_rt, 0.1_rt};
        pp.addarr("tke_values", tke_vals);
        pp.add("surface_temp_flux", 0.0_rt);
    }
    // Transport
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }

    // Initialize necessary parts of solver
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Create turbulence model
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    // Get turbulence model
    auto& tmodel = sim().turbulence_model();

    // Get coefficients
    auto model_dict = tmodel.model_coeffs();

    // Constants for fields
    const amrex::Real srate = 20.0_rt;
    const amrex::Real Tgz = 0.0_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real x3 = 1016.0_rt;
    const amrex::Real lscale_s = (lambda * kappa * x3) / (lambda + kappa * x3);
    const amrex::Real tlscale_val = lscale_s;
    const amrex::Real tke_val = 0.1_rt;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_strain_field(vel, srate);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    // Set up temperature field with constant gradient in z
    auto& temp = sim().repo().get_field("temperature");
    init_temperature_field(temp, Tgz);
    // Give values to tlscale and tke arrays
    auto& tlscale = sim().repo().get_field("turb_lscale");
    tlscale.setVal(tlscale_val);
    auto& tke = sim().repo().get_field("tke");
    tke.setVal(tke_val);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    const auto max_val = utils::field_max(muturb);
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real epsilon = kynema_sgf::utils::powi(Cmu, 3) *
                                std::pow(tke_val, 1.5_rt) /
                                (tlscale_val + 1.0e-3_rt);
    const amrex::Real stratification = 0.0_rt;
    const amrex::Real Rt =
        kynema_sgf::utils::powi(tke_val / epsilon, 2) * stratification;
    const amrex::Real Cmu_Rt =
        (0.556_rt + 0.108_rt * Rt) /
        (1.0_rt + 0.308_rt * Rt + 0.00837_rt * kynema_sgf::utils::powi(Rt, 2));
    const amrex::Real tol = 0.12_rt;
    const amrex::Real nut_max =
        rho0 * Cmu_Rt * tlscale_val * std::sqrt(tke_val);
    EXPECT_NEAR(max_val, nut_max, tol);
}

// ForestDrag.canopy_tke adds fd (beta_p |U|^3 - beta_d |U| k) to the
// KransAxell source in the canopy only; off, a forest changes nothing
TEST_F(TurbRANSTest, test_1eqKrans_canopy_tke_source)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().init_physics();
    sim().create_transport_model();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& repo = sim().repo();
    auto& vel = repo.get_field("velocity");
    init_strain_field(vel, 0.01_rt);
    repo.get_field("density").setVal(rho0);
    init_temperature_field(repo.get_field("temperature"), 0.0_rt);
    repo.get_field("turb_lscale").setVal(10.0_rt);
    const amrex::Real tke_val = 0.1_rt;
    repo.get_field("tke").setVal(tke_val);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const amrex::Real dt = 0.5_rt;
    sim().time().delta_t() = dt;

    auto& src_no_forest = repo.declare_field("src_no_forest", 1, 0, 1);
    auto& src_off = repo.declare_field("src_off", 1, 0, 1);
    auto& src_on = repo.declare_field("src_on", 1, 0, 1);
    src_no_forest.setVal(0.0_rt);
    src_off.setVal(0.0_rt);
    src_on.setVal(0.0_rt);

    // Default inputs, before the forest exists
    kynema_sgf::pde::tke::KransAxell krans_off(sim());
    krans_off(0, kynema_sgf::FieldState::New, src_no_forest(0));

    // Canopy below 100 m
    const amrex::Real fd = 0.02_rt;
    auto& forest = repo.declare_field("forest_drag", 1, 1, 1);
    set_below(forest, fd, 100.0_rt);
    krans_off(0, kynema_sgf::FieldState::New, src_off(0));

    {
        amrex::ParmParse pp("ForestDrag");
        pp.add("canopy_tke", true);
    }
    kynema_sgf::pde::tke::KransAxell krans_on(sim());
    krans_on(0, kynema_sgf::FieldState::New, src_on(0));

    // Off: bit-identical to no forest
    amrex::MultiFab::Subtract(src_off(0), src_no_forest(0), 0, 0, 1, 0);
    EXPECT_EQ(src_off(0).norm0(), 0.0_rt);

    // On: the canopy term in the canopy (z = 88 m), nothing above (z = 104 m)
    amrex::MultiFab::Subtract(src_on(0), src_no_forest(0), 0, 0, 1, 0);
    const int i = 5;
    const int j = 5;
    const int kc = 5;
    const amrex::Real ux = utils::field_probe(vel, 0, i, j, kc, 0);
    const amrex::Real uy = utils::field_probe(vel, 0, i, j, kc, 1);
    const amrex::Real uz = utils::field_probe(vel, 0, i, j, kc, 2);
    const amrex::Real ws = std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
    const amrex::Real expected =
        (fd * ws * ws * ws) +
        (tke_val * std::expm1(-4.0_rt * fd * ws * dt) / dt);
    EXPECT_NEAR(
        utils::field_probe(src_on, 0, i, j, kc), expected,
        1.0e-6_rt * std::abs(expected));
    EXPECT_EQ(utils::field_probe(src_on, 0, i, j, kc + 1), 0.0_rt);
    EXPECT_EQ(utils::field_probe(src_on, 0, i, j, 40), 0.0_rt);
}

// KLAxell with the forest present and ForestDrag.canopy_tke off: the eddy
// viscosity is bit-identical to the run without a forest
TEST_F(TurbRANSTest, test_1eqKrans_canopy_length_off)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().init_physics();
    sim().create_transport_model();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& repo = sim().repo();
    init_strain_field(repo.get_field("velocity"), 0.01_rt);
    repo.get_field("density").setVal(rho0);
    init_temperature_field(repo.get_field("temperature"), 0.0_rt);
    repo.get_field("turb_lscale").setVal(10.0_rt);
    repo.get_field("tke").setVal(0.1_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto& mu_no_forest = repo.declare_field("mu_no_forest", 1, 0, 1);
    amrex::MultiFab::Copy(
        mu_no_forest(0), repo.get_field("mu_turb")(0), 0, 0, 1, 0);

    set_below(repo.declare_field("forest_drag", 1, 1, 1), 0.02_rt, 100.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    amrex::MultiFab::Subtract(
        mu_no_forest(0), repo.get_field("mu_turb")(0), 0, 0, 1, 0);
    EXPECT_EQ(mu_no_forest(0).norm0(), 0.0_rt);
}

// ForestDrag.canopy_tke on: l = min(l_s, alpha / fd) in the canopy only
TEST_F(TurbRANSTest, test_1eqKrans_canopy_length_on)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    {
        amrex::ParmParse pp("ForestDrag");
        pp.add("canopy_tke", true);
    }
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().init_physics();
    sim().create_transport_model();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& repo = sim().repo();
    init_strain_field(repo.get_field("velocity"), 0.01_rt);
    repo.get_field("density").setVal(rho0);
    init_temperature_field(repo.get_field("temperature"), 0.0_rt);
    repo.get_field("turb_lscale").setVal(10.0_rt);
    const amrex::Real tke_val = 0.1_rt;
    repo.get_field("tke").setVal(tke_val);

    // Zero forest drag: no limit anywhere
    auto& forest = repo.declare_field("forest_drag", 1, 1, 1);
    forest.setVal(0.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto& mu_free = repo.declare_field("mu_free", 1, 0, 1);
    amrex::MultiFab::Copy(mu_free(0), repo.get_field("mu_turb")(0), 0, 0, 1, 0);

    // Canopy below 100 m: alpha / fd = 2 m, below l_s at every height
    const amrex::Real fd = 0.02_rt;
    set_below(forest, fd, 100.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);

    const auto& mu = repo.get_field("mu_turb");
    const auto& tlscale = repo.get_field("turb_lscale");
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real lcap = 0.04_rt / fd;
    const amrex::Real mu_cap = rho0 * Cmu * lcap * std::sqrt(tke_val);
    const amrex::Real tol = 10.0_rt * kynema_sgf::constants::TIGHT_TOL;
    for (int k = 0; k < 6; ++k) {
        EXPECT_NEAR(utils::field_probe(tlscale, 0, 5, 5, k), lcap, tol);
        EXPECT_NEAR(utils::field_probe(mu, 0, 5, 5, k), mu_cap, tol);
    }
    // Above the canopy (z = 104 m and up) the eddy viscosity is unchanged
    for (int k = 6; k < 64; k += 7) {
        EXPECT_EQ(
            utils::field_probe(mu, 0, 5, 5, k),
            utils::field_probe(mu_free, 0, 5, 5, k));
    }
}

namespace {
// Relative tolerance for values computed along two paths (device vs host
// pow, l + EPS vs l): tight in double, a few ulp in single precision
constexpr amrex::Real rtol =
    std::is_same_v<amrex::Real, float> ? 1.0e-5_rt : 1.0e-10_rt;

// Set up KLAxell + KransAxell after the inputs are in ParmParse. When
// zlo_wall is given, the domain spans [zlo_wall, zlo_wall + 1024] in z.
void setup_klaxell(
    kynema_sgf::CFDSim& sim,
    const amrex::Real rho0,
    const amrex::Real tke_val,
    const amrex::Real tgrad)
{
    sim.pde_manager().register_icns();
    sim.init_physics();
    sim.create_transport_model();
    sim.create_turbulence_model();
    sim.turbulence_model().post_init_actions();
    auto& repo = sim.repo();
    init_strain_field(repo.get_field("velocity"), 0.01_rt);
    repo.get_field("density").setVal(rho0);
    init_temperature_field(repo.get_field("temperature"), tgrad);
    repo.get_field("turb_lscale").setVal(10.0_rt);
    repo.get_field("tke").setVal(tke_val);
}

amrex::Real neutral_lscale(const amrex::Real kappa, const amrex::Real z)
{
    const amrex::Real lambda = 30.0_rt;
    return (lambda * kappa * z) / (lambda + kappa * z);
}
} // namespace

// KransAxell builds dissipation from KLAxell_coeffs.Cmu, and KLAxell builds
// its length scale from ABL.kappa. Before, KransAxell used ABL.Cmu (default
// 0.556) and KLAxell a fixed kappa = 0.41.
TEST_F(TurbRANSTest, test_1eqKrans_shared_cmu_kappa)
{
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real Cmu = 0.5_rt;
    const amrex::Real kappa = 0.4_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    {
        amrex::ParmParse pp("KLAxell_coeffs");
        pp.add("Cmu", Cmu);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("kappa", kappa);
    }
    populate_parameters();
    initialize_mesh();
    const amrex::Real tke_val = 0.1_rt;
    setup_klaxell(sim(), rho0, tke_val, 0.0_rt);
    auto& repo = sim().repo();
    sim().turbulence_model().update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    sim().time().delta_t() = 0.5_rt;
    auto& src = repo.declare_field("src_cmu", 1, 0, 1);
    src.setVal(0.0_rt);
    kynema_sgf::pde::tke::KransAxell krans(sim());
    krans(0, kynema_sgf::FieldState::New, src(0));

    const auto& tlscale = repo.get_field("turb_lscale");
    const auto& dissip = repo.get_field("dissipation");
    const amrex::Real dz = 16.0_rt;
    for (int k = 1; k < 64; k += 9) {
        const amrex::Real l = neutral_lscale(kappa, (k + 0.5_rt) * dz);
        EXPECT_NEAR(utils::field_probe(tlscale, 0, 5, 5, k), l, rtol * l);
        const amrex::Real eps = Cmu * Cmu * Cmu * std::pow(tke_val, 1.5_rt) /
                                (l + kynema_sgf::constants::EPS);
        EXPECT_NEAR(utils::field_probe(dissip, 0, 5, 5, k), eps, rtol * eps);
    }
}

// The flat-ground length scale and the wall value of k are measured from
// prob_lo[2], not from z = 0
TEST_F(TurbRANSTest, test_1eqKrans_wall_distance_problo)
{
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real zlo = 100.0_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, zlo}};
        amrex::Vector<amrex::Real> probhi{
            {1024.0_rt, 1024.0_rt, zlo + 1024.0_rt}};
        pp.addarr("prob_lo", problo);
        pp.addarr("prob_hi", probhi);
    }
    initialize_mesh();
    ASSERT_EQ(mesh().Geom(0).ProbLo(2), zlo);
    const amrex::Real tke_val = 0.1_rt;
    setup_klaxell(sim(), rho0, tke_val, 0.0_rt);
    auto& repo = sim().repo();
    sim().turbulence_model().update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const amrex::Real dt = 0.5_rt;
    sim().time().delta_t() = dt;
    auto& src = repo.declare_field("src_problo", 1, 0, 1);
    src.setVal(0.0_rt);
    kynema_sgf::pde::tke::KransAxell krans(sim());
    krans(0, kynema_sgf::FieldState::New, src(0));

    const amrex::Real kappa = 0.41_rt;
    const amrex::Real dz = 16.0_rt;
    const auto& tlscale = repo.get_field("turb_lscale");
    for (int k = 0; k < 4; ++k) {
        const amrex::Real l = neutral_lscale(kappa, (k + 0.5_rt) * dz);
        EXPECT_NEAR(utils::field_probe(tlscale, 0, 5, 5, k), l, rtol * l);
    }

    // First cell: production - dissipation + (k_wall - k) / (5 dt), with
    // k_wall = u*^2 / Cmu^2 and u* from the cell above at 3 dz / 2
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real z0 = 0.1_rt;
    const auto& vel = repo.get_field("velocity");
    const amrex::Real ux = utils::field_probe(vel, 0, 5, 5, 1, 0);
    const amrex::Real uy = utils::field_probe(vel, 0, 5, 5, 1, 1);
    const amrex::Real ustar =
        std::sqrt((ux * ux) + (uy * uy)) * kappa / std::log(1.5_rt * dz / z0);
    const amrex::Real k_wall = ustar * ustar / (Cmu * Cmu);
    const amrex::Real expected =
        utils::field_probe(repo.get_field("shear_prod"), 0, 5, 5, 0) +
        utils::field_probe(repo.get_field("buoy_prod"), 0, 5, 5, 0) -
        utils::field_probe(repo.get_field("dissipation"), 0, 5, 5, 0) +
        ((k_wall - tke_val) / (5.0_rt * dt));
    EXPECT_NEAR(
        utils::field_probe(src, 0, 5, 5, 0), expected,
        rtol * std::abs(expected));
}

namespace {
// k = 0 in a stratified run: Rt was 0/0 and mu_t NaN. Checks mu_t, the
// production terms and the thermal diffusivity from update_alphaeff.
void check_zero_tke_stratified(kynema_sgf::CFDSim& sim, const amrex::Real rho0)
{
    // Stable gradient, so N^2 > 0 and the Rt = (k / eps)^2 N^2 branch runs
    setup_klaxell(sim, rho0, 0.0_rt, 0.01_rt);
    sim.turbulence_model().update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto& repo = sim.repo();
    EXPECT_FALSE(repo.get_field("mu_turb")(0).contains_nan());
    EXPECT_FALSE(repo.get_field("buoy_prod")(0).contains_nan());
    EXPECT_FALSE(repo.get_field("shear_prod")(0).contains_nan());
    EXPECT_FALSE(repo.get_field("turb_lscale")(0).contains_nan());
    EXPECT_EQ(utils::field_max(repo.get_field("mu_turb")), 0.0_rt);
    // The thermal diffusivity from update_alphaeff, with l > 0 as left by an
    // earlier step (k = 0 and l = 0 together is also caught by epsilon = 0/0)
    repo.get_field("turb_lscale").setVal(10.0_rt);
    auto& alpha = repo.declare_cc_field("alpha_test", 1, 1, 1);
    alpha.set_default_fillpatch_bc(sim.time());
    alpha.setVal(0.0_rt);
    sim.turbulence_model().update_alphaeff(alpha);
    EXPECT_FALSE(alpha(0).contains_nan());
}

// Non-neutral, so the neutral switch does not set Rt = 0
void add_nonneutral_flux()
{
    amrex::ParmParse pp("ABL");
    pp.add("surface_temp_flux", -0.05_rt);
}
} // namespace

TEST_F(TurbRANSTest, test_1eqKrans_zero_tke_stratified)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    add_nonneutral_flux();
    populate_parameters();
    initialize_mesh();
    check_zero_tke_stratified(sim(), rho0);
}

// Same with KLAxell_coeffs.legacy_alphaeff = true, whose update_alphaeff
// recomputes Rt and carries its own k = 0 guard
TEST_F(TurbRANSTest, test_1eqKrans_zero_tke_stratified_legacy)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    add_nonneutral_flux();
    {
        amrex::ParmParse pp("KLAxell_coeffs");
        pp.add("legacy_alphaeff", true);
    }
    populate_parameters();
    initialize_mesh();
    check_zero_tke_stratified(sim(), rho0);
}

// ABL.horizontal_sponge_tke on flat ground: k is set to the profile value in
// one step inside the east sponge and unchanged outside. Before, the sponge
// ran only with TerrainDrag.
TEST_F(TurbRANSTest, test_1eqKrans_flat_horizontal_sponge)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    initialize_mesh();
    const amrex::Real tke_val = 0.3_rt;
    setup_klaxell(sim(), rho0, tke_val, 0.0_rt);
    auto& repo = sim().repo();
    sim().turbulence_model().update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const amrex::Real dt = 0.5_rt;
    sim().time().delta_t() = dt;

    auto& src_off = repo.declare_field("src_sponge_off", 1, 0, 1);
    auto& src_on = repo.declare_field("src_sponge_on", 1, 0, 1);
    src_off.setVal(0.0_rt);
    src_on.setVal(0.0_rt);
    kynema_sgf::pde::tke::KransAxell krans_off(sim());
    krans_off(0, kynema_sgf::FieldState::New, src_off(0));
    {
        amrex::ParmParse pp("ABL");
        pp.add("horizontal_sponge_tke", true);
    }
    {
        amrex::ParmParse pp("DragForcing");
        pp.add("sponge_east", true);
        pp.add("sponge_distance_east", 200.0_rt);
    }
    kynema_sgf::pde::tke::KransAxell krans_on(sim());
    krans_on(0, kynema_sgf::FieldState::New, src_on(0));
    amrex::MultiFab::Subtract(src_on(0), src_off(0), 0, 0, 1, 0);

    // Profile k = 0.1 at every height; the sponge starts at x = 824 m. At
    // and below meso_sponge_start (600 m, k <= 36) only the horizontal
    // sponge acts; above it only the vertical sponge, so the difference is 0.
    const amrex::Real expected = -(tke_val - 0.1_rt) / dt;
    for (int k = 40; k < 64; k += 11) {
        for (int i = 0; i < 10; ++i) {
            EXPECT_EQ(utils::field_probe(src_on, 0, i, 5, k), 0.0_rt);
        }
    }
    for (int k = 2; k < 37; k += 11) {
        for (int i = 8; i < 10; ++i) {
            EXPECT_NEAR(
                utils::field_probe(src_on, 0, i, 5, k), expected,
                rtol * std::abs(expected));
        }
        for (int i = 0; i < 8; ++i) {
            EXPECT_EQ(utils::field_probe(src_on, 0, i, 5, k), 0.0_rt);
        }
    }
}

// ABL.meso_sponge_start defaults to 600 m in KLAxell, as in KransAxell and
// the docs (KLAxell used 2000 m). Neutral run, stable gradient: l = l_s at
// and below 600 m, the stratified length scale above.
TEST_F(TurbRANSTest, test_1eqKrans_meso_start_default)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    initialize_mesh();
    setup_klaxell(sim(), rho0, 0.1_rt, 0.01_rt);
    sim().turbulence_model().update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& tlscale = sim().repo().get_field("turb_lscale");
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real dz = 16.0_rt;
    // z = 584 m and 648 m
    const amrex::Real l_below = neutral_lscale(kappa, (36 + 0.5_rt) * dz);
    const amrex::Real l_above = neutral_lscale(kappa, (40 + 0.5_rt) * dz);
    EXPECT_NEAR(
        utils::field_probe(tlscale, 0, 5, 5, 36), l_below, rtol * l_below);
    EXPECT_LT(utils::field_probe(tlscale, 0, 5, 5, 40), 0.9_rt * l_above);
}

namespace {
// Neutral run with a stable gradient; returns alpha_eff - alpha_lam and
// mu_t at height index k
std::pair<amrex::Real, amrex::Real>
klaxell_alpha_t(kynema_sgf::CFDSim& sim, const amrex::Real rho0, const int k)
{
    setup_klaxell(sim, rho0, 0.1_rt, 0.01_rt);
    sim.turbulence_model().update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto& repo = sim.repo();
    auto& alpha = repo.declare_cc_field("alpha_t_test", 1, 1, 1);
    alpha.set_default_fillpatch_bc(sim.time());
    alpha.setVal(0.0_rt);
    sim.turbulence_model().update_alphaeff(alpha);
    const auto lam = sim.transport_model().alpha();
    // Constant laminar diffusivity
    const amrex::Real alpha_lam = (*lam)(0).max(0);
    EXPECT_EQ((*lam)(0).min(0), alpha_lam);
    const amrex::Real alpha_t =
        utils::field_probe(alpha, 0, 5, 5, k) - alpha_lam;
    return {alpha_t, utils::field_probe(repo.get_field("mu_turb"), 0, 5, 5, k)};
}
} // namespace

// alpha_t = mu_t / Pr_t(Rt) with the Rt of update_turbulent_viscosity. Below
// meso_sponge_start in a neutral run Rt = 0, so alpha_t = mu_t, matching the
// buoyancy production built from C_mu'(0) = C_mu.
TEST_F(TurbRANSTest, test_1eqKrans_alphaeff_same_rt)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    initialize_mesh();
    const auto [alpha_t, mu_t] = klaxell_alpha_t(sim(), rho0, 10);
    ASSERT_GT(mu_t, 0.0_rt);
    EXPECT_NEAR(alpha_t, mu_t, rtol * mu_t);
}

// Above meso_sponge_start (z = 648 m) Rt > 0, and alpha_t = mu_t / Pr_t(Rt)
// with the Rt of the viscosity kernel. That Rt uses the length scale on entry
// (10 m, set by setup_klaxell).
TEST_F(TurbRANSTest, test_1eqKrans_alphaeff_stratified_rt)
{
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real Tref = 265.0_rt;
    const amrex::Real gravz = 10.0_rt;
    add_klaxell_inputs(Tref, gravz, rho0);
    populate_parameters();
    initialize_mesh();
    const auto [alpha_t, mu_t] = klaxell_alpha_t(sim(), rho0, 40);
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real tke = 0.1_rt;
    const amrex::Real l_in = 10.0_rt;
    const amrex::Real N2 = gravz / Tref * 0.01_rt;
    const amrex::Real eps = Cmu * Cmu * Cmu * std::pow(tke, 1.5_rt) / l_in;
    const amrex::Real Rt = (tke / eps) * (tke / eps) * N2;
    const amrex::Real Pr = (1.0_rt + 0.193_rt * Rt) / (1.0_rt + 0.0302_rt * Rt);
    ASSERT_GT(mu_t, 0.0_rt);
    ASSERT_GT(Pr, 1.5_rt);
    EXPECT_NEAR(alpha_t, mu_t / Pr, rtol * mu_t);
}

// KLAxell_coeffs.legacy_alphaeff = true keeps the earlier form: Rt
// recomputed from the current l, without the neutral switch
TEST_F(TurbRANSTest, test_1eqKrans_alphaeff_legacy)
{
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real Tref = 265.0_rt;
    const amrex::Real gravz = 10.0_rt;
    add_klaxell_inputs(Tref, gravz, rho0);
    {
        amrex::ParmParse pp("KLAxell_coeffs");
        pp.add("legacy_alphaeff", true);
    }
    populate_parameters();
    initialize_mesh();
    const int k = 10;
    const auto [alpha_t, mu_t] = klaxell_alpha_t(sim(), rho0, k);
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real tke = 0.1_rt;
    const amrex::Real l = neutral_lscale(0.41_rt, (k + 0.5_rt) * 16.0_rt);
    const amrex::Real N2 = gravz / Tref * 0.01_rt;
    const amrex::Real eps = Cmu * Cmu * Cmu * std::pow(tke, 1.5_rt) / l;
    const amrex::Real Rt = (tke / eps) * (tke / eps) * N2;
    ASSERT_GT(Rt, 0.0_rt);
    const amrex::Real Pr = (1.0_rt + 0.193_rt * Rt) / (1.0_rt + 0.0302_rt * Rt);
    EXPECT_NEAR(alpha_t, mu_t / Pr, rtol * mu_t);
}

TEST_F(TurbRANSTest, test_1eqKrans_production_per_unit_mass)
{
    // Stable column with surface_temp_flux = 0 and every cell below
    // ABL.meso_sponge_start (set to 2000 m): KLAxell uses l = l_s(z) and Rt =
    // 0, so Cmu(Rt) = Cmu'(Rt) = Cmu and, per unit mass,
    //   shear_prod = S^2 nu_t,  buoy_prod = -nu_t N^2,  nu_t = mu_t / rho
    const amrex::Real Tref = 300.0_rt;
    const amrex::Real gravz = 9.81_rt;
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real srate = 2.0_rt;
    const amrex::Real Tgz = 0.003_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "KLAxell");
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -gravz};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_rate", 0.0_rt);
        // Every cell (z < 1024 m) below the neutral switch height; the
        // KLAxell default is 600 m
        pp.add("meso_sponge_start", 2000.0_rt);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 4000.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{300.0_rt, 312.0_rt};
        pp.addarr("temperature_values", t_vals);
        pp.add("surface_temp_flux", 0.0_rt);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }

    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& vel = sim().repo().get_field("velocity");
    init_strain_field(vel, srate);
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    auto& temp = sim().repo().get_field("temperature");
    init_temperature_field(temp, Tgz);
    auto& tlscale = sim().repo().get_field("turb_lscale");
    auto& tke = sim().repo().get_field("tke");

    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real tke_val = 0.5_rt;
    const amrex::Real nsq = gravz * Tgz / Tref;

    // Largest |computed / expected - 1| over the domain for mu_t and both
    // production terms
    auto max_errors = [&]() {
        const auto& shear_prod = sim().repo().get_field("shear_prod");
        const auto& buoy_prod = sim().repo().get_field("buoy_prod");
        const auto& muturb = sim().repo().get_field("mu_turb");
        const auto& problo = sim().mesh().Geom(0).ProbLoArray();
        const amrex::Real dz = sim().mesh().Geom(0).CellSize()[2];
        const auto& sp_arrs = shear_prod(0).const_arrays();
        const auto& bp_arrs = buoy_prod(0).const_arrays();
        const auto& mu_arrs = muturb(0).const_arrays();
        const auto r = amrex::ParReduce(
            amrex::TypeList<
                amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real, amrex::Real, amrex::Real>{},
            shear_prod(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real> {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dz);
                const amrex::Real ls =
                    (lambda * kappa * z) / (lambda + (kappa * z));
                const amrex::Real nut = Cmu * ls * std::sqrt(tke_val);
                const amrex::Real mu_ref = rho0 * nut;
                const amrex::Real sp_ref = srate * srate * nut;
                const amrex::Real bp_ref = -nut * nsq;
                return {
                    std::abs((mu_arrs[nbx](i, j, k) - mu_ref) / mu_ref),
                    std::abs((sp_arrs[nbx](i, j, k) - sp_ref) / sp_ref),
                    std::abs((bp_arrs[nbx](i, j, k) - bp_ref) / bp_ref)};
            });
        std::array<amrex::Real, 3> out{
            amrex::get<0>(r), amrex::get<1>(r), amrex::get<2>(r)};
        amrex::ParallelDescriptor::ReduceRealMax(out.data(), 3);
        return out;
    };

    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e3_rt;
    tke.setVal(tke_val);
    tlscale.setVal(1.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto err = max_errors();
    EXPECT_NEAR(err[0], 0.0_rt, tol);
    EXPECT_NEAR(err[1], 0.0_rt, tol);
    EXPECT_NEAR(err[2], 0.0_rt, tol);

    // Same check through the terrain kernel: flat terrain at z = 0 with no
    // blanked cells gives the same z, l_s, mu_t and production
    {
        auto& blank = sim().repo().declare_int_field("terrain_blank", 1, 1, 1);
        blank.setVal(0);
        auto& height = sim().repo().declare_field("terrain_height", 1, 1, 1);
        height.setVal(0.0_rt);
    }
    tke.setVal(tke_val);
    tlscale.setVal(1.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    err = max_errors();
    EXPECT_NEAR(err[0], 0.0_rt, tol);
    EXPECT_NEAR(err[1], 0.0_rt, tol);
    EXPECT_NEAR(err[2], 0.0_rt, tol);
}

} // namespace kynema_sgf_tests
