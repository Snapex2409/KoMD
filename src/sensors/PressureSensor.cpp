//
// Created by alex on 8/8/24.
//

#include "PressureSensor.h"
#include "Registry.h"
#include "util/constants.h"
#include "potentials/Limit.h"
#include "potentials/ATM2B.h"
#include "potentials/FENE.h"
#include "potentials/LJ12_6.h"
#include "potentials/ATM_NOLIST.h"
#include "potentials/ATM.h"
#include "IO/Logging.h"
#include <type_traits>
#include <utility>
#include <sstream>
#include <fstream>


PressureSensor::PressureSensor() : Sensor("PressureSensor"),
m_nu(Registry::instance->configuration()->energy_3b)
{
    m_pressures = KW::vec_t<double>("Pressure Tensor", NUM_PRESSURES);
    m_nu = m_nu * 1e+9 * Constants::conv_J_Ei; // convert external unit to internal
}

void PressureSensor::measure()
{
    for (int idx = 0; idx < NUM_PRESSURES; ++idx) m_pressures[idx] = 0;
    contributeMomentum();
    contributePotentials();
}

void PressureSensor::write(uint64_t simstep)
{
    std::stringstream file_name;
    file_name << "pressures_" << simstep << ".txt";
    std::ofstream file(file_name.str());
    if (!file.is_open()) {
        Log::io->error() << "Could not write file: " << file_name.str() << std::endl;
        return;
    }
    // what we actually have stored atm is the stress tensor -> need to div by volume
    // v * v * m = (m² * kg) / s²
    // force unit is u * A / ps²
    // pressure unit is u / (ps² * A)
    // output unit will be Pa = N/m² = kg * m/s² * m^-2 = kg / (s² * m)
    const SciValue convA_m {1.0, -10};
    const SciValue convps_s {1.0, -12};
    const SciValue factor = Constants::conv_Da_kg / (convps_s * convps_s * convA_m);
    const double volume = (Registry::instance->configuration()->domainHigh - Registry::instance->configuration()->domainLow).product();

    int next_count = 1;
    int count = 0;
    for (int idx = 0; idx < NUM_PRESSURES; ++idx)
    {
        file << m_pressures[idx] / volume * factor;
        count++;

        if (count == next_count)
        {
            count = 0;
            next_count = next_count + 1;
            file << std::endl;
        }
        else file << " ";
    }
    file.flush();
    file.close();
}

void PressureSensor::contributeMomentum()
{
    auto container = Registry::instance->moleculeContainer();
    auto size = container->getNumMolecules();

    std::array<double, NUM_PRESSURES> c{0};
    Kokkos::parallel_reduce("Pressure Sensor - Momentum", size,
        MomentumKernel(container->getSOA().v(), container->getSOA().mass()),
        c[0], c[1], c[2], c[3], c[4], c[5]);
    Kokkos::fence("Pressure Sensor - Momentum");
    for (int idx = 0; idx < NUM_PRESSURES; ++idx) m_pressures[idx] += c[idx];
}

void PressureSensor::contributePotentials()
{
    auto& potentials = Registry::instance->forceFunctors();
    auto& potentials3B = Registry::instance->forceFunctors3b();
    for (auto ff_ptr : potentials) runFF2(ff_ptr);
    for (auto ff_ptr : potentials3B) runFF3(ff_ptr);
}

void PressureSensor::runFF2(std::shared_ptr<ForceFunctor> ff)
{
    if (auto ptr = std::dynamic_pointer_cast<Limit>(ff); ptr != nullptr) return;
    if (auto ptr = std::dynamic_pointer_cast<FENE>(ff); ptr != nullptr) return;

    auto container = Registry::instance->moleculeContainer();
    auto& soa = container->getSOA();
    auto& pairList = container->getPairList();
    const double cutoff2 = std::pow(Registry::instance->configuration()->cutoff, 2);

    std::array<double, NUM_PRESSURES> c{0};
    if (auto ptr = std::dynamic_pointer_cast<LJ12_6>(ff); ptr != nullptr)
    {
        Kokkos::parallel_reduce("Pressure Sensor - Potential 2B", pairList.size(),
        LJ12_6_Pressure(pairList.getPairs(), pairList.getOffsets(), soa.f(), soa.id(), soa.r(), soa.sigma(), soa.epsilon(), cutoff2),
        c[0], c[1], c[2], c[3], c[4], c[5]);
    }
    if (auto ptr = std::dynamic_pointer_cast<ATM2B>(ff); ptr != nullptr)
    {
        Kokkos::parallel_reduce("Pressure Sensor - Potential 2B", pairList.size(),
        ATM2B_Pressure(pairList.getPairs(), pairList.getOffsets(), soa.f(), soa.id(), soa.r(), soa.sigma(), soa.epsilon(), cutoff2, ptr->getPotentialFactor()),
        c[0], c[1], c[2], c[3], c[4], c[5]);
    }

    Kokkos::fence("Pressure Sensor - Potential 2B");
    for (int idx = 0; idx < NUM_PRESSURES; ++idx) m_pressures[idx] += c[idx];
}

void PressureSensor::runFF3(std::shared_ptr<ForceFunctor3B> ff)
{
    if (auto ptr = std::dynamic_pointer_cast<ATM_NOLIST>(ff); ptr != nullptr) return;

    auto container = Registry::instance->moleculeContainer();
    auto& soa = container->getSOA();
    auto& tripleList = container->getTripleList();
    const double cutoff2 = std::pow(Registry::instance->configuration()->cutoff, 2);

    std::array<double, NUM_PRESSURES> c{0};
    if (auto ptr = std::dynamic_pointer_cast<ATM>(ff); ptr != nullptr)
    {
        Kokkos::parallel_reduce("Pressure Sensor - Potential 3B", tripleList.size(),
        ATM_Pressure(tripleList.getTriplets(), tripleList.getOffsets(), soa.f(), soa.r(), m_nu, cutoff2),
        c[0], c[1], c[2], c[3], c[4], c[5]);
    }

    Kokkos::fence("Pressure Sensor - Potential 3B");
    for (int idx = 0; idx < NUM_PRESSURES; ++idx) m_pressures[idx] += c[idx];
}

void PressureSensor::MomentumKernel::operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy, double& acc_xz, double& acc_yz, double& acc_zz) const
{
    acc_xx += v[idx][0] * v[idx][0] * m[idx];
    acc_xy += v[idx][0] * v[idx][1] * m[idx];
    acc_yy += v[idx][1] * v[idx][1] * m[idx];
    acc_xz += v[idx][0] * v[idx][2] * m[idx];
    acc_yz += v[idx][1] * v[idx][2] * m[idx];
    acc_zz += v[idx][2] * v[idx][2] * m[idx];
}

void PressureSensor::LJ12_6_Pressure::operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy,
                                                 double& acc_xz, double& acc_yz, double& acc_zz) const
{
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);

    if (id(s_idx_0) == id(s_idx_1)) return;
    const math::d3 shift0 = pair_offsets(idx, 0);
    const math::d3 shift1 = pair_offsets(idx, 1);

    const math::d3 dr = (r[s_idx_0] + shift0) - (r[s_idx_1] + shift1);
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps[s_idx_0] * eps[s_idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const auto force = dr * fac;

    //if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], force);
    //if (shift1 == 0) Kokkos::atomic_sub(&f[s_idx_1], force);
    // for now assuming we use all forces
    acc_xx += 0.5 * (dr.x() * force.x() + dr.x() * force.x());
    acc_xy += 0.5 * (dr.x() * force.y() + dr.y() * force.x());
    acc_yy += 0.5 * (dr.y() * force.y() + dr.y() * force.y());
    acc_xz += 0.5 * (dr.x() * force.z() + dr.z() * force.x());
    acc_yz += 0.5 * (dr.y() * force.z() + dr.z() * force.y());
    acc_zz += 0.5 * (dr.z() * force.z() + dr.z() * force.z());
}

void PressureSensor::ATM2B_Pressure::operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy, double& acc_xz,
                                                double& acc_yz, double& acc_zz) const
{
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);

    if (id(s_idx_0) == id(s_idx_1)) return;
    const math::d3 shift0 = pair_offsets(idx, 0);
    const math::d3 shift1 = pair_offsets(idx, 1);

    const math::d3 dr = (r[s_idx_0] + shift0) - (r[s_idx_1] + shift1);
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps[s_idx_0] * eps[s_idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const auto force = dr * fac * potential_factor;

    //if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], force);
    //if (shift1 == 0) Kokkos::atomic_sub(&f[s_idx_1], force);
    // for now assuming we use all forces
    acc_xx += 0.5 * (dr.x() * force.x() + dr.x() * force.x());
    acc_xy += 0.5 * (dr.x() * force.y() + dr.y() * force.x());
    acc_yy += 0.5 * (dr.y() * force.y() + dr.y() * force.y());
    acc_xz += 0.5 * (dr.x() * force.z() + dr.z() * force.x());
    acc_yz += 0.5 * (dr.y() * force.z() + dr.z() * force.y());
    acc_zz += 0.5 * (dr.z() * force.z() + dr.z() * force.z());
}

void PressureSensor::ATM_Pressure::operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy, double& acc_xz,
                                              double& acc_yz, double& acc_zz) const
{
    const uint64_t s_idx_0 = triplets(idx, 0);
    const uint64_t s_idx_1 = triplets(idx, 1);
    const uint64_t s_idx_2 = triplets(idx, 2);

    const math::d3 shift0 = triple_offsets(idx, 0);
    const math::d3 shift1 = triple_offsets(idx, 1);
    const math::d3 shift2 = triple_offsets(idx, 2);

    const math::d3 r0 = r[s_idx_0] + shift0;
    const math::d3 r1 = r[s_idx_1] + shift1;
    const math::d3 r2 = r[s_idx_2] + shift2;

    const math::d3 dr0_1 = r0 - r1;
    const math::d3 dr0_2 = r0 - r2;
    const math::d3 dr1_2 = r1 - r2;

    const double dr0_1_2 = dr0_1.dot(dr0_1);
    const double dr0_2_2 = dr0_2.dot(dr0_2);
    const double dr1_2_2 = dr1_2.dot(dr1_2);
    // symmetrical condition
    if (dr0_1_2 > cutoff2 || dr0_2_2 > cutoff2 || dr1_2_2 > cutoff2) return;

    // we are within range
    // compute forces
    const double r0_1 = Kokkos::sqrt(dr0_1_2);
    const double r0_2 = Kokkos::sqrt(dr0_2_2);
    const double r1_2 = Kokkos::sqrt(dr1_2_2);

    const double dVdR01 = - ATM::ATM_Force::comp_force(nu, r0_1, r0_2, r1_2);
    const double dVdR02 = - ATM::ATM_Force::comp_force(nu, r0_2, r0_1, r1_2);
    const double dVdR12 = - ATM::ATM_Force::comp_force(nu, r1_2, r0_1, r0_2);
    const math::d3 f01 = dr0_1 * dVdR01;
    const math::d3 f02 = dr0_2 * dVdR02;
    const math::d3 f12 = dr1_2 * dVdR12;
    // for now assuming we use all forces
    // ijk === 012
    acc_xx += 0.5 * (
        dr0_1.x() * f01.x() + dr0_1.x() * f01.x() +
        dr1_2.x() * f12.x() + dr1_2.x() * f12.x() +
        dr0_2.x() * f02.x() + dr0_2.x() * f02.x()
    );
    acc_xy += 0.5 * (
        dr0_1.x() * f01.y() + dr0_1.y() * f01.x() +
        dr1_2.x() * f12.y() + dr1_2.y() * f12.x() +
        dr0_2.x() * f02.y() + dr0_2.y() * f02.x()
    );
    acc_yy += 0.5 * (
        dr0_1.y() * f01.y() + dr0_1.y() * f01.y() +
        dr1_2.y() * f12.y() + dr1_2.y() * f12.y() +
        dr0_2.y() * f02.y() + dr0_2.y() * f02.y()
    );
    acc_xz += 0.5 * (
        dr0_1.x() * f01.z() + dr0_1.z() * f01.x() +
        dr1_2.x() * f12.z() + dr1_2.z() * f12.x() +
        dr0_2.x() * f02.z() + dr0_2.z() * f02.x()
    );
    acc_yz += 0.5 * (
        dr0_1.y() * f01.z() + dr0_1.z() * f01.y() +
        dr1_2.y() * f12.z() + dr1_2.z() * f12.y() +
        dr0_2.y() * f02.z() + dr0_2.z() * f02.y()
    );
}