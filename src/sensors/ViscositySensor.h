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
        TensorEvaluator(bool origin_free, double dt, std::initializer_list<PressureSensor::PressureDim> eval_points);
        void eval(const KW::vec_t<double>& tensor);
        double correlate(int origin_begin, int origin_end, int origin_stride, int offset, int buffer_idx, double function_shift = 0.0) const;
        double integrate() const;
        double ens_avg_pressure(int buffer_idx) const;
    private:
        int m_num_tensor_points = 0;
        std::array<bool, PressureSensor::PressureDim::NUM_PRESSURES> m_tensor_points {};

        std::vector<std::vector<double>> m_pressures;
        int m_num_evals = 0;

        int m_correlation_stride = 1;
        double m_dt = 0;
        bool m_origin_free = false;
        double m_volume = 0;
    };

    TensorEvaluator m_bulk_evaluator;
    TensorEvaluator m_shear_evaluator;
};



#endif //VISCOSITYSENSOR_H
