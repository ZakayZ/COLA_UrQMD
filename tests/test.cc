#include <Urqmd4/Urqmd4Module.hh>
#include <COLA/EventData.hh>

#include "expected_particles_data.hh"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>

namespace {

void expectRelNear(double a, double b, double relTol) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    EXPECT_NEAR(a, b, relTol * scale);
}

} // namespace

TEST(UrQmdParticleTest, MatchesColaReference) {
    auto mod = cola::fortran::Urqmd4Module();
    auto filters = mod.GetModuleFilters();
    std::unordered_map<std::string, std::string> params = {
        {"config_file", "urqmd_input"}
    };
    auto filter = filters["URQMDGenerator"]->Create(params);
    auto* gen = dynamic_cast<cola::VGenerator*>(filter.get());
    ASSERT_NE(gen, nullptr);

    auto data = (*gen)();
    ASSERT_NE(data, nullptr);

    const auto& parts = data->particles;
    ASSERT_EQ(parts.size(), kExpectedParticleCount);

    constexpr double kTol = 1e-7;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        SCOPED_TRACE("particle index: " + std::to_string(i));
        const auto& p = parts[i];
        const auto& e = kExpectedParticles[i];
        const auto& pos = p.position;
        const auto& mom = p.momentum;

        expectRelNear(pos.t, e.t, kTol);
        expectRelNear(pos.x, e.x, kTol);
        expectRelNear(pos.y, e.y, kTol);
        expectRelNear(pos.z, e.z, kTol);
        expectRelNear(mom.e, e.e, kTol);
        expectRelNear(mom.x, e.px, kTol);
        expectRelNear(mom.y, e.py, kTol);
        expectRelNear(mom.z, e.pz, kTol);
        EXPECT_EQ(p.pdgCode, e.pdg);
        EXPECT_EQ(static_cast<int>(p.pClass), e.pclass);
    }
}
