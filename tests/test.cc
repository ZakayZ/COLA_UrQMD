#include "pdg_quantum_numbers.hh"

#include <COLA.hh>
#include <COLA/EventData.hh>
#include <COLA_UrQMD/COLA_UrQMDModule.hh>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>

namespace {

  constexpr double kRelTol = 1e-7;
  constexpr int kCollisionsPerEnergy = 4;

  double Tolerance(double a, double b, double rel_tolerance = kRelTol) {
    return rel_tolerance * std::max({1.0, std::abs(a), std::abs(b)});
  }

  std::unique_ptr<cola::VGenerator> MakeUrQMDGenerator(double elb_a_gev, int random_seed,
                                                       const std::string& generated_config_path) {
    auto module = cola::fortran::COLA_UrQMDModule();
    auto filters = module.GetModuleFilters();
    std::ostringstream elb_stream;
    elb_stream << elb_a_gev;

    std::unordered_map<std::string, std::string> params = {
        {"pro", "197 79"},
        {"tar", "197 79"},
        {"nev", "1"},
        {"imp", "5."},
        {"elb", elb_stream.str()},
        {"tim", "200 200"},
        {"rsd", std::to_string(random_seed)},
        {"tables_file", "tables.dat"},
        {"generated_config_file", generated_config_path},
    };
    auto filter = filters["URQMDGenerator"]->Create(params);

    return std::unique_ptr<cola::VGenerator>(dynamic_cast<cola::VGenerator*>(filter.release()));
  }

}  // namespace

class UrQmdConservationTest : public ::testing::TestWithParam<double> {};

TEST_P(UrQmdConservationTest, MomentumBaryonAndChargeConservation) {
  const double lab_energy_per_nucleon = GetParam();
  const auto config_path =
      std::string("urqmd_conservation_") + std::to_string(static_cast<int>(lab_energy_per_nucleon)) + ".txt";
  auto generator =
      MakeUrQMDGenerator(lab_energy_per_nucleon, 91000 + static_cast<int>(lab_energy_per_nucleon), config_path);
  ASSERT_NE(generator, nullptr);

  for (int collision_index = 0; collision_index < kCollisionsPerEnergy; ++collision_index) {
    SCOPED_TRACE("collision: " + std::to_string(collision_index) +
                 ", elb AGeV: " + std::to_string(lab_energy_per_nucleon));
    auto event_data = (*generator)();
    ASSERT_NE(event_data, nullptr);

    const auto initial_four_momentum = std::accumulate(
        event_data->ini_state.ini_state_particles.begin(), event_data->ini_state.ini_state_particles.end(),
        cola::LorentzVector{}, [](cola::LorentzVector acc, const cola::Particle& particle) {
          acc += particle.momentum;
          return acc;
        });
    const auto final_four_momentum =
        std::accumulate(event_data->particles.begin(), event_data->particles.end(), cola::LorentzVector{},
                        [](cola::LorentzVector acc, const cola::Particle& particle) {
                          acc += particle.momentum;
                          return acc;
                        });

    EXPECT_NEAR(initial_four_momentum.e, final_four_momentum.e,
                Tolerance(initial_four_momentum.e, final_four_momentum.e, 5e-4))
        << "energy mismatch (lab frame)";
    EXPECT_NEAR(initial_four_momentum.x, final_four_momentum.x,
                Tolerance(initial_four_momentum.x, final_four_momentum.x, 5e-4))
        << "px mismatch (lab frame)";
    EXPECT_NEAR(initial_four_momentum.y, final_four_momentum.y,
                Tolerance(initial_four_momentum.y, final_four_momentum.y, 5e-4))
        << "py mismatch (lab frame)";
    EXPECT_NEAR(initial_four_momentum.z, final_four_momentum.z,
                Tolerance(initial_four_momentum.z, final_four_momentum.z, 5e-4))
        << "pz mismatch (lab frame)";

    const uint64_t initial_baryon_number = std::accumulate(
        event_data->ini_state.ini_state_particles.begin(), event_data->ini_state.ini_state_particles.end(), 0LL,
        [](uint64_t sum, const cola::Particle& particle) { return sum + urqmd_test::BaryonNumber(particle.pdg_code); });
    const uint64_t final_baryon_number = std::accumulate(
        event_data->particles.begin(), event_data->particles.end(), 0LL,
        [](uint64_t sum, const cola::Particle& particle) { return sum + urqmd_test::BaryonNumber(particle.pdg_code); });
    EXPECT_EQ(initial_baryon_number, final_baryon_number);

    const uint64_t initial_charge_times3 = std::accumulate(
        event_data->ini_state.ini_state_particles.begin(), event_data->ini_state.ini_state_particles.end(), 0LL,
        [](uint64_t sum, const cola::Particle& particle) { return sum + urqmd_test::HepChgTimes3(particle.pdg_code); });
    const uint64_t final_charge_times3 = std::accumulate(
        event_data->particles.begin(), event_data->particles.end(), 0LL,
        [](uint64_t sum, const cola::Particle& particle) { return sum + urqmd_test::HepChgTimes3(particle.pdg_code); });
    EXPECT_EQ(initial_charge_times3, final_charge_times3);
  }
}

INSTANTIATE_TEST_SUITE_P(LabEnergyPerNucleon, UrQmdConservationTest, ::testing::Values(15.0, 40.0, 100.0, 160.0));
