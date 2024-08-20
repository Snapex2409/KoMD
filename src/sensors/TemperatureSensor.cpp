//
// Created by alex on 8/8/24.
//

#include "TemperatureSensor.h"
#include "Registry.h"
#include "math/SciValue.h"
#include "math/Geometry.h"
#include "util/constants.h"

TemperatureSensor::TemperatureSensor() : Sensor("Temperature"),
p_low(Registry::instance->configuration()->domainLow),
p_high(Registry::instance->configuration()->domainHigh),
m_enable_opt_approx(Registry::instance->configuration()->sensor_temp_opt),
m_temperature(0) { }

TemperatureSensor::TemperatureSensor(const math::d3 &low, const math::d3 &high) :
Sensor("Temperature"), p_low(low), p_high(high),
m_enable_opt_approx(Registry::instance->configuration()->sensor_temp_opt),
m_temperature(0) { }

void TemperatureSensor::measure() {
    auto container = Registry::instance->moleculeContainer();

    SOA::vec_t<double> mv2("MV2", 1);
    SOA::vec_scatter_t<double> mv2_scatter(mv2);
    SOA::vec_t<double> num_sites("Num_Sites", 1);
    SOA::vec_scatter_t<double> num_sites_scatter(num_sites);
    mv2[0] = 0;
    num_sites[0] = 0;

    for (auto it = container->iteratorCell(MoleculeContainer::DOMAIN); it.isValid(); ++it) {
        auto& cell = it.cell();
        auto& soa = cell.soa();
        Kokkos::parallel_for("Temp Measurement", soa.size(), Temperature_Kernel(soa.r(), soa.v(), soa.mass(), p_low, p_high, mv2_scatter, num_sites_scatter));
    }
    Kokkos::fence("Temp fence");
    Kokkos::Experimental::contribute(mv2, mv2_scatter);
    Kokkos::Experimental::contribute(num_sites, num_sites_scatter);

    static const SciValue convDaInvkb = Constants::conv_Da_kg / Constants::kB;
    static const SciValue factor = Constants::conv_Aps_ms * Constants::conv_Aps_ms * convDaInvkb;
    num_sites[0] *= 3; // for 3 dimensions
    m_temperature = factor * (mv2[0] / num_sites[0]);
}

void TemperatureSensor::write(uint64_t simstep) { }

double TemperatureSensor::getTemperature() { return m_temperature; }

void TemperatureSensor::Temperature_Kernel::operator()(int idx) const {
    if (!math::pointInBox(r[idx], low, high)) return;
    const auto vel = v[idx];

    auto mv2_access = mv2.access();
    auto num_sites_access = num_sites.access();

    mv2_access(0) += m[idx] * vel.dot(vel);
    num_sites_access(0) += 1;
}
