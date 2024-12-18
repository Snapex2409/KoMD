//
// Created by alex on 16.12.24.
//

#include "ViscositySensor.h"

#include "Registry.h"
#include "util/constants.h"
#include "IO/Logging.h"
#include <utility>
#include <sstream>
#include <fstream>

using PSD = PressureSensor::PressureDim;

ViscositySensor::ViscositySensor() : Sensor("Viscosity"), m_pressure_sensor(),
                                     m_sample_window_size(Registry::instance->configuration()->sensor_visc_window),
m_bulk_evaluator(Registry::instance->configuration()->sensor_visc_no_orig,
                 Registry::instance->configuration()->delta_t,
                 {PSD::XX, PSD::YY, PSD::ZZ}),
m_shear_evaluator(Registry::instance->configuration()->sensor_visc_no_orig,
                  Registry::instance->configuration()->delta_t,
                 {PSD::XY, PSD::XZ, PSD::YZ}),
m_dt(Registry::instance->configuration()->delta_t)
{
}

void ViscositySensor::measure()
{
    m_pressure_sensor.measure();
    auto& pressures = m_pressure_sensor.getPressureTensor();
    m_bulk_evaluator.eval(pressures);
    m_shear_evaluator.eval(pressures);
}

void ViscositySensor::write(uint64_t simstep)
{
    // Units:
    // Pa: kg / (s² * m)
    // square: kg² / (s^4 * m²)
    // int-t:  kg² / (s^3 * m²)
    // rescale:kg / (s * m)
    // require: 2x conv to Pa 1x conv to s

    // V/(kB*T) * int_t(Pa*Pa)
    // m³/J * int_t(kg²/(s^4 * m²))
    // m³/(N*m) * kg²/(s^3 * m²)
    // m²/(kg * m / s²) * kg²/(s^3 * m²)
    // m*s²/kg * kg²/(s^3 * m²)
    // kg/(s * m)

    // we have:
    // V/(kB*T) * int_t(P*P)
    // A³/J * int_t(IP*IP) =
    // A³/J * int_t(u²/(ps^4 * A²)) =
    // A³/J * u²/(ps^3 * A²) =
    // (A*u²)/(ps³*J)

    // stress: (kg * m²) / s²
    // square: (kg² * m^4) / s^4
    // int-t:  (kg² * m^4) / s^3
    // rescale:(kg * m²) / s
    const SciValue convA_m {1.0, -10};
    const SciValue convps_s {1.0, -12};
    const SciValue conv = (Constants::conv_Da_kg * Constants::conv_Da_kg * convA_m) / (convps_s * convps_s * convps_s);
    const double V = (Registry::instance->configuration()->domainHigh - Registry::instance->configuration()->domainLow).product();

    const SciValue factor = conv * (V / (Constants::kB * Registry::instance->configuration()->temperature));

    std::stringstream file_name;
    file_name << "viscosity_" << simstep << ".txt";
    std::ofstream file(file_name.str());
    if (!file.is_open()) {
        Log::io->error() << "Could not write file: " << file_name.str() << std::endl;
        return;
    }

    file << "b: " << m_bulk_evaluator.integrate() * factor << std::endl;
    file << "s: " << m_shear_evaluator.integrate() * factor << std::endl;
    file.flush();
    file.close();
}

ViscositySensor::TensorEvaluator::TensorEvaluator(bool origin_free, double dt,
    std::initializer_list<PressureSensor::PressureDim> eval_points) : m_num_evals(0), m_correlation_stride(1),
    m_dt(dt), m_origin_free(origin_free)
{
    m_pressures.resize(PressureSensor::PressureDim::NUM_PRESSURES);

    m_num_tensor_points = eval_points.size();
    for (const auto p : eval_points) m_tensor_points[p] = true;

    m_volume = (Registry::instance->configuration()->domainHigh - Registry::instance->configuration()->domainLow).product();
}

void ViscositySensor::TensorEvaluator::eval(const KW::vec_t<double>& tensor)
{
    for (int idx = 0; idx < m_tensor_points.size(); ++idx)
    {
        if (m_tensor_points[idx]) m_pressures[idx].push_back(tensor[idx] / m_volume);
    }
    m_num_evals++;
}

double ViscositySensor::TensorEvaluator::correlate(int origin_begin, int origin_end, int origin_stride, int offset,
    int buffer_idx, double function_shift) const
{
    const int num_steps = (origin_end - origin_begin) / origin_stride;
    double tmp = 0.0;
    for (int origin = origin_begin; origin < origin_end; origin += origin_stride)
    {
        tmp += (m_pressures[buffer_idx][origin] - function_shift) * (m_pressures[buffer_idx][origin + offset] - function_shift);
    }
    tmp /= num_steps;
    return tmp;
}

double ViscositySensor::TensorEvaluator::integrate() const
{
    std::vector<double> integrals;

    for (int pDim = 0; pDim < PressureSensor::PressureDim::NUM_PRESSURES; ++pDim)
    {
        if (!m_tensor_points[pDim]) continue;

        double shift = 0;
        if (m_origin_free) shift = ens_avg_pressure(pDim);

        // first compute function values for all integration nodes
        std::vector<double> correlations(m_num_evals, 0.0);
        for (int idx = 0; idx < m_num_evals; idx++)
        {
            correlations[idx] = correlate(0, m_num_evals - idx, m_correlation_stride, idx, pDim, shift);
        }

        // using TS here
        double integral = 0;
        integral += correlations[0] + correlations[m_num_evals-1];
        for (int idx = 1; idx < m_num_evals-1; idx++) integral += correlations[idx];
        integral *= m_dt / 2.0;

        integrals.push_back(integral);
    }

    // page 274
    double sum = 0;
    for (int idx = 0; idx < integrals.size(); idx++) sum += integrals[idx];
    return sum / integrals.size();
}

double ViscositySensor::TensorEvaluator::ens_avg_pressure(int buffer_idx) const
{
    double avg = 0;
    for (int idx = 0; idx < m_num_evals; idx++)
    {
        avg += m_pressures[buffer_idx][idx];
    }
    return avg / m_num_evals;
}
