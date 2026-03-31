#include <Urqmd4/Urqmd4Module.hh>
#include <COLA/EventData.hh>

#include "expected_particles_data.hh"
#include "pdg_quantum_numbers.hh"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>

namespace {

constexpr double RelTol = 1e-7;
constexpr int CollisionsPerEnergy = 4;
constexpr double GeVToMeV = 1000.0;

double Tolerance(double a, double b, double relTolerance = RelTol) {
    return relTolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

std::unique_ptr<cola::VGenerator> MakeUrQMDGenerator(double elbAGeV, int randomSeed,
                                                     const std::string& generatedConfigPath) {
    auto module = cola::fortran::Urqmd4Module();
    auto filters = module.GetModuleFilters();
    std::ostringstream elbStream;
    elbStream << elbAGeV;

    std::unordered_map<std::string, std::string> params = {
        {"pro", "197 79"},
        {"tar", "197 79"},
        {"nev", "1"},
        {"imp", "5."},
        {"elb", elbStream.str()},
        {"tim", "200 200"},
        {"rsd", std::to_string(randomSeed)},
        {"tables_file", "tables.dat"},
        {"generated_config_file", generatedConfigPath},
    };
    auto filter = filters["URQMDGenerator"]->Create(params);

    return std::unique_ptr<cola::VGenerator>(dynamic_cast<cola::VGenerator*>(filter.release()));
}

} // namespace

TEST(UrQmdParticleTest, MatchesColaReference) {
    auto module = cola::fortran::Urqmd4Module();
    auto filters = module.GetModuleFilters();
    std::unordered_map<std::string, std::string> params = {
        {"config_file", "urqmd_input"}
    };
    auto filter = filters["URQMDGenerator"]->Create(params);
    auto* generator = dynamic_cast<cola::VGenerator*>(filter.get());
    ASSERT_NE(generator, nullptr);

    auto eventData = (*generator)();
    ASSERT_NE(eventData, nullptr);

    const auto& particles = eventData->particles;
    ASSERT_EQ(particles.size(), ExpectedParticleCount);

    for (std::size_t index = 0; index < particles.size(); ++index) {
        SCOPED_TRACE("particle index: " + std::to_string(index));
        const auto& particle = particles[index];
        const auto& expected = ExpectedParticles[index];
        const auto expectedMomentum = expected.GetMomentum() * GeVToMeV;
        const auto expectedPosition = expected.GetPosition();
        const auto& position = particle.position;
        const auto& momentum = particle.momentum;

        EXPECT_NEAR(position.t, expectedPosition.t, Tolerance(position.t, expectedPosition.t, RelTol));
        EXPECT_NEAR(position.x, expectedPosition.x, Tolerance(position.x, expectedPosition.x, RelTol));
        EXPECT_NEAR(position.y, expectedPosition.y, Tolerance(position.y, expectedPosition.y, RelTol));
        EXPECT_NEAR(position.z, expectedPosition.z, Tolerance(position.z, expectedPosition.z, RelTol));
        EXPECT_NEAR(momentum.e, expectedMomentum.e, Tolerance(momentum.e, expectedMomentum.e, RelTol));
        EXPECT_NEAR(momentum.x, expectedMomentum.x, Tolerance(momentum.x, expectedMomentum.x, RelTol));
        EXPECT_NEAR(momentum.y, expectedMomentum.y, Tolerance(momentum.y, expectedMomentum.y, RelTol));
        EXPECT_NEAR(momentum.z, expectedMomentum.z, Tolerance(momentum.z, expectedMomentum.z, RelTol));
        EXPECT_EQ(particle.pdgCode, expected.pdg);
        EXPECT_EQ(static_cast<int>(particle.pClass), expected.pclass);
    }
}

class UrQmdConservationTest : public ::testing::TestWithParam<double> {};

TEST_P(UrQmdConservationTest, MomentumBaryonAndChargeConservation) {
    const double labEnergyPerNucleon = GetParam();
    const auto configPath =
        std::string("urqmd_conservation_") + std::to_string(static_cast<int>(labEnergyPerNucleon)) + ".txt";
    auto generator = MakeUrQMDGenerator(labEnergyPerNucleon, 91000 + static_cast<int>(labEnergyPerNucleon), configPath);
    ASSERT_NE(generator, nullptr);

    for (int collisionIndex = 0; collisionIndex < CollisionsPerEnergy; ++collisionIndex) {
        SCOPED_TRACE("collision: " + std::to_string(collisionIndex) +
                     ", elb AGeV: " + std::to_string(labEnergyPerNucleon));
        auto eventData = (*generator)();
        ASSERT_NE(eventData, nullptr);

        const auto initialFourMomentum =
            std::accumulate(eventData->iniState.iniStateParticles.begin(),
                            eventData->iniState.iniStateParticles.end(), cola::LorentzVector{},
                            [](cola::LorentzVector acc, const cola::Particle& particle) {
                                acc += particle.momentum;
                                return acc;
                            });
        const auto finalFourMomentum =
            std::accumulate(eventData->particles.begin(), eventData->particles.end(), cola::LorentzVector{},
                            [](cola::LorentzVector acc, const cola::Particle& particle) {
                                acc += particle.momentum;
                                return acc;
                            });

        EXPECT_NEAR(initialFourMomentum.e, finalFourMomentum.e,
                    Tolerance(initialFourMomentum.e, finalFourMomentum.e, 5e-4))
            << "energy mismatch (lab frame)";
        EXPECT_NEAR(initialFourMomentum.x, finalFourMomentum.x,
                    Tolerance(initialFourMomentum.x, finalFourMomentum.x, 5e-4))
            << "px mismatch (lab frame)";
        EXPECT_NEAR(initialFourMomentum.y, finalFourMomentum.y,
                    Tolerance(initialFourMomentum.y, finalFourMomentum.y, 5e-4))
            << "py mismatch (lab frame)";
        EXPECT_NEAR(initialFourMomentum.z, finalFourMomentum.z,
                    Tolerance(initialFourMomentum.z, finalFourMomentum.z, 5e-4))
            << "pz mismatch (lab frame)";

        const uint64_t initialBaryonNumber =
            std::accumulate(eventData->iniState.iniStateParticles.begin(),
                            eventData->iniState.iniStateParticles.end(), 0LL,
                            [](uint64_t sum, const cola::Particle& particle) {
                                return sum + urqmd_test::BaryonNumber(particle.pdgCode);
                            });
        const uint64_t finalBaryonNumber =
            std::accumulate(eventData->particles.begin(), eventData->particles.end(), 0LL,
                            [](uint64_t sum, const cola::Particle& particle) {
                                return sum + urqmd_test::BaryonNumber(particle.pdgCode);
                            });
        EXPECT_EQ(initialBaryonNumber, finalBaryonNumber);

        const uint64_t initialChargeTimes3 =
            std::accumulate(eventData->iniState.iniStateParticles.begin(),
                            eventData->iniState.iniStateParticles.end(), 0LL,
                            [](uint64_t sum, const cola::Particle& particle) {
                                return sum + urqmd_test::HepChgTimes3(particle.pdgCode);
                            });
        const uint64_t finalChargeTimes3 =
            std::accumulate(eventData->particles.begin(), eventData->particles.end(), 0LL,
                            [](uint64_t sum, const cola::Particle& particle) {
                                return sum + urqmd_test::HepChgTimes3(particle.pdgCode);
                            });
        EXPECT_EQ(initialChargeTimes3, finalChargeTimes3);
    }
}

INSTANTIATE_TEST_SUITE_P(
    LabEnergyPerNucleon,
    UrQmdConservationTest,
    ::testing::Values(15.0, 40.0, 100.0, 160.0));
