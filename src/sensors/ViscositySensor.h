//
// Created by alex on 16.12.24.
//

#ifndef VISCOSITYSENSOR_H
#define VISCOSITYSENSOR_H

#include "PressureSensor.h"
#include "Sensor.h"

class ViscositySensor : public Sensor {
public:
    ViscositySensor();
    ~ViscositySensor() override = default;
    void measure() override;
    void write(uint64_t simstep) override;
    /// @returns results converted back to standard SI units
    double getBulkViscosity();
    /// @returns results converted back to standard SI units
    double getShearViscosity();

private:
    PressureSensor m_pressure_sensor;
    /// steps per ens avg
    const int m_sample_window_size;
    /// delta t
    const double m_dt;

    class TensorEvaluator
    {
    public:
        TensorEvaluator(int window_size, bool origin_free, double dt, std::initializer_list<PressureSensor::PressureDim> eval_points);
        void eval(const KW::vec_t<double>& tensor);
        double integrate() const;
        void reset() { m_ens_averages.clear(); m_ens_averages_times.clear(); }

    private:
        bool m_origin_free = false;
        bool m_origin_init = false;
        double m_origin = 0;

        int m_window_size = 0;
        int m_sample_counter = 0;

        int m_num_tensor_points = 0;
        std::array<bool, PressureSensor::PressureDim::NUM_PRESSURES> m_tensor_points {};

        double m_averaging_buffer;
        std::vector<double> m_ens_averages;
        std::vector<double> m_ens_averages_times;
        double m_zero_quantity = 0;
        double m_dt = 0;
    };

    TensorEvaluator m_bulk_evaluator;
    TensorEvaluator m_shear_evaluator;
};



#endif //VISCOSITYSENSOR_H
