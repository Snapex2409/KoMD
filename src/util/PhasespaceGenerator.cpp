//
// Created by alex on 8/4/24.
//

#include "PhasespaceGenerator.h"
#include "util/constants.h"
#include "Registry.h"
#include "molecule/Molecule.h"

#include <random>

math::d3 random_rotate(double alpha, double beta, double gamma, const math::d3& vec) {
    double x = std::cos(beta) * std::cos(gamma) * vec.x() +
               (std::sin(alpha) * std::sin(beta) * std::cos(gamma) - std::cos(alpha) * std::sin(gamma)) * vec.y() +
               (std::cos(alpha) * std::sin(beta) * std::cos(gamma) + std::sin(alpha) * std::sin(gamma)) * vec.z();
    double y = std::cos(beta) * std::sin(gamma) * vec.x() +
               (std::sin(alpha) * std::sin(beta) * std::sin(gamma) + std::cos(alpha) * std::cos(gamma)) * vec.y() +
               (std::cos(alpha) * std::sin(beta) * std::sin(gamma) - std::sin(alpha) * std::cos(gamma)) * vec.z();
    double z = -std::sin(beta) * vec.x() +
               std::sin(alpha) * std::cos(beta) * vec.y() +
               std::cos(alpha) * std::cos(beta) * vec.z();
    return {x, y, z};
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-flp30-c"
void PhasespaceGenerator::generate() {
    // for now we just create a fixed tetrahedral molecule with fixed epsilon and sigma
    double density = Registry::instance->configuration()->density;
    auto container = Registry::instance->moleculeContainer();
    auto config = Registry::instance->configuration();

    static const std::array<math::d3, 4> coords {
        math::d3 {0.35355339, 0.35355339, 0.35355339},
        math::d3 {-0.35355339, -0.35355339, 0.35355339},
        math::d3 {-0.35355339, 0.35355339, -0.35355339},
        math::d3 {0.35355339, -0.35355339, -0.35355339}
    };

    uint64_t counter = 0;
    const auto spacing = static_cast<uint64_t>(1.0 / density);
    static std::default_random_engine rng(43); // NOLINT(*-msc51-cpp) I know... that's the point
    std::uniform_real_distribution<double> uniform(0, M_PI * 2);

    for (double z = config->domainLow.z() + 0.5; z < config->domainHigh.z(); z += 1.0) {
        for (double y = config->domainLow.y() + 0.5; y < config->domainHigh.y(); y += 1.0) {
            for (double x = config->domainLow.x() + 0.5; x < config->domainHigh.x(); x += 1.0) {
                if (counter % spacing == 0) counter++;
                else {
                    counter++;
                    continue;
                }

                const math::d3 center_pos {x, y, z};
                Molecule molecule;
                double alpha = uniform(rng), beta = uniform(rng), gamma = uniform(rng);

                for (int idx = 0; idx < 4; idx++) {
                    const math::d3 site_pos = random_rotate(alpha, beta, gamma, center_pos + coords[idx]);
                    const math::d3 v = getMaxwellBoltzmannVelocity(config->temperature, 1.0);
                    molecule.addSite(config->epsilon, config->sigma, 1.0, site_pos, v);
                }
            }
        }
    }
}
#pragma clang diagnostic pop

math::d3 PhasespaceGenerator::getMaxwellBoltzmannVelocity(double T, double m) {
    static std::default_random_engine rng(42); // NOLINT(*-msc51-cpp) I know... that's the point
    static const SciValue kbInvConvDa = Constants::kB / Constants::conv_Da_kg;

    double std_dev = std::sqrt((T / m) * kbInvConvDa);
    std::normal_distribution<double> normal(0, 1);
    math::d3 result {normal(rng), normal(rng), normal(rng)};
    result *= std_dev;
    return result;
}
