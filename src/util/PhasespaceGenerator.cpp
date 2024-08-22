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
    auto& ps_regions = config->phasespace_gen_regions;
    auto& components = Registry::instance->components();

    const math::d3 domain_size = config->domainHigh - config->domainLow;
    const double domain_volume = domain_size.product();
    const double num_molecules = density * domain_volume;
    const double mol_per_dim = std::pow(num_molecules, 1./3.);
    const math::d3 spacing = domain_size / mol_per_dim;

    static std::default_random_engine rng(43); // NOLINT(*-msc51-cpp) I know... that's the point
    std::uniform_real_distribution<double> uniform(0, M_PI * 2);

    for (auto [begin, end, cid] : ps_regions) {
        for (double z = begin.z() + 0.5; z < end.z(); z += spacing.z()) {
            for (double y = begin.y() + 0.5; y < end.y(); y += spacing.y()) {
                for (double x = begin.x() + 0.5; x < end.x(); x += spacing.x()) {
                    const math::d3 center_pos {x, y, z};
                    Molecule molecule;
                    molecule.setCID(cid);
                    Component& component = components[cid];

                    //velocity and rotation should be same for all sites
                    double alpha = uniform(rng), beta = uniform(rng), gamma = uniform(rng);
                    const math::d3 v = getMaxwellBoltzmannVelocity(config->temperature, 1.0);
                    for (int s_idx = 0; s_idx < component.getSites().size(); s_idx++) {
                        Site& site = component.getSites()[s_idx];
                        const math::d3 site_pos = center_pos + random_rotate(alpha, beta, gamma, site.r_arr());
                        molecule.addSite(site.getEpsilon(), site.getSigma(), site.getMass(), site_pos, v);
                    }

                    container->addMolecule(molecule);
                }
            }
        }
    }
}
#pragma clang diagnostic pop

math::d3 PhasespaceGenerator::getMaxwellBoltzmannVelocity(double T, double m) {
    static std::default_random_engine rng(42); // NOLINT(*-msc51-cpp) I know... that's the point
    static const SciValue kbInvConvDa = Constants::kB / Constants::conv_Da_kg;

    double std_dev = std::sqrt((T / m) * kbInvConvDa); // this is in m/s
    std_dev *= Constants::conv_ms_Aps; // now in A/ps (base unit of simulation)
    std::normal_distribution<double> normal(0, 1);
    math::d3 result {normal(rng), normal(rng), normal(rng)};
    result *= std_dev;
    return result;
}
