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
m_bulk_evaluator(Registry::instance->configuration()->sensor_visc_window,
                 Registry::instance->configuration()->sensor_visc_no_orig,
                 Registry::instance->configuration()->delta_t,
                 {PSD::XX, PSD::YY, PSD::ZZ}),
m_shear_evaluator(Registry::instance->configuration()->sensor_visc_window,
                 Registry::instance->configuration()->sensor_visc_no_orig,
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
    const SciValue convA_m {1.0, -10};
    const SciValue convps_s {1.0, -12};
    const SciValue convIP_Pa = Constants::conv_Da_kg / (convps_s * convps_s * convA_m);
    const SciValue factor = convIP_Pa * convIP_Pa / convps_s;

    std::stringstream file_name;
    file_name << "viscosity_" << simstep << ".txt";
    std::ofstream file(file_name.str());
    if (!file.is_open()) {
        Log::io->error() << "Could not write file: " << file_name.str() << std::endl;
        return;
    }

    file << "b: " << m_bulk_evaluator.integrate() * factor << std::endl;
    file << "s: " << m_shear_evaluator.integrate() * factor << std::endl;
    //m_bulk_evaluator.reset();
    //m_shear_evaluator.reset();
    file.flush();
    file.close();
}

ViscositySensor::TensorEvaluator::TensorEvaluator(int window_size, bool origin_free, double dt,
    std::initializer_list<PressureSensor::PressureDim> eval_points) : m_window_size(window_size), m_origin_free(origin_free), m_dt(dt)
{
    m_origin_init = false;
    m_origin = 0;

    m_num_tensor_points = eval_points.size();
    for (const auto p : eval_points) m_tensor_points[p] = true;

    m_averaging_buffer = 0;
    m_zero_quantity = 0;
}

void ViscositySensor::TensorEvaluator::eval(const KW::vec_t<double>& tensor)
{
    m_sample_counter++;

    double tmp = 0;
    for (int pos = 0; pos < m_tensor_points.size(); pos++)
    {
        if (m_tensor_points[pos]) tmp += tensor[pos];
    };
    tmp /= m_num_tensor_points;

    // init origin free
    if (m_origin_free && !m_origin_init)
    {
        m_averaging_buffer += tmp;

        if (m_sample_counter % m_window_size == 0)
        {
            m_origin = m_averaging_buffer / m_window_size;
            m_origin_init = true;
            m_averaging_buffer = 0;
            m_zero_quantity = tmp - m_origin;
            m_sample_counter = 0;
        }
        return;
    }

    // init with origin
    if (!m_origin_free && !m_origin_init)
    {
        m_origin_init = true;
        m_averaging_buffer = 0;
        m_zero_quantity = tmp;
        m_sample_counter = 0;
        return;
    }

    // handle main measurement origin free
    if (m_origin_free)
    {
        m_averaging_buffer += (tmp - m_origin) * m_zero_quantity;

        if (m_sample_counter % m_window_size == 0)
        {
            m_ens_averages.push_back(m_averaging_buffer / m_window_size);
            m_ens_averages_times.push_back(Registry::instance->simulation()->simstep() * m_dt);
            m_sample_counter = 0;
            m_averaging_buffer = 0;
        }
    }
    // handle main measurement with origin
    else
    {
        m_averaging_buffer += tmp * m_zero_quantity;

        if (m_sample_counter % m_window_size == 0)
        {
            m_ens_averages.push_back(m_averaging_buffer / m_window_size);
            m_ens_averages_times.push_back(Registry::instance->simulation()->simstep() * m_dt);
            m_sample_counter = 0;
            m_averaging_buffer = 0;
        }
    }
}

double ViscositySensor::TensorEvaluator::integrate() const
{
    // using TS here
    double integral = 0;
    const double delta_t = m_window_size * m_dt;
    integral += m_ens_averages[0] + m_ens_averages[m_ens_averages.size()-1];
    for (int idx = 1; idx < m_ens_averages.size()-1; idx++) integral += m_ens_averages[idx];
    integral *= delta_t / 2.0;

    return integral;
}
