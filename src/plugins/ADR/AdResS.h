//
// Created by alex on 10/15/24.
//

#ifndef ADRESS_H
#define ADRESS_H

#include "plugins/Plugin.h"
#include "math/Array.h"
#include "util/Kokkos_Wrapper.h"

class AdResS : public Plugin {
public:
    AdResS();

    ~AdResS() override = default;

    void init() override;

    void pre_main_loop() override;

    void pre_container_update() override;

    void post_container_update() override;

    void post_forces() override;

    /// computes the weight for the given position
    KOKKOS_FUNCTION static double weight(const math::d3 &r, const math::d3 &fp_low, const math::d3 &fp_high, const math::d3 &h_dim);

    /// Kernel for device
    struct Weight_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> coms;
        KW::vec_t<double> weights;

        math::d3 fp_low;
        math::d3 fp_high;
        math::d3 h_dim;
    };

    /// Kernel for device
    struct Force_Distribute_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> f;
        KW::vec_t<double> m;
        const int site_count;
        const double total_mass;
    };

    struct CG_Pos_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> r;
        KW::vec_t<math::d3> com;
        KW::vec_t<double> m;
    };

private:
    /// FP region low
    const math::d3 m_fp_low;
    /// FP region high
    const math::d3 m_fp_high;
    /// H region dimensions (excluding FP)
    const math::d3 m_h_dim;
    /// buffer for weights, populated during force comp, used in integration
    KW::vec_t<double> m_weights;
    /// hybrid total mass
    double m_total_mass;
    /// total site count per molecule
    int m_site_count;
    /// flag to check if main loop is running
    bool m_is_first;
};



#endif //ADRESS_H
