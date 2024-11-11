//
// Created by alex on 11/5/24.
//

#ifndef IBI_FORCEFUNCTOR_H
#define IBI_FORCEFUNCTOR_H

#include "IBI_Math.h"
#include "potentials/ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"
#include <type_traits>

struct IBI_Reload {};
struct IBI_Default {};

template<typename T, bool reload = std::is_same_v<IBI_Reload, T>>
requires(std::is_same_v<T, IBI_Reload> or std::is_same_v<T, IBI_Default>)
class IBI_ForceFunctor : public ForceFunctor {
public:
    IBI_ForceFunctor();
    ~IBI_ForceFunctor() override = default;

    FunctionPL& getForceFunction() { return m_force; }
    FunctionPL& getPotentialFunction() { return m_potential; }

    /// Kernel for device
    struct IBI_Force {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of PairList::pairs
        KW::nvec_t<int, 2> pairs;
        /// shallow copy of PairList::pair_offsets
        KW::nvec_t<math::d3, 2> pair_offsets;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of FunctionPL::x_values
        KW::vec_t<double> force_x_values;
        /// shallow copy of FunctionPL::y_values
        KW::vec_t<double> force_y_values;
        /// default lower bound value of FunctionPL
        double force_def_low;
        /// default uppoer bound value of FunctionPL
        double force_def_high;
        /// cutoff radius squared
        const double cutoff2;
        /// exclusion zone for reload: low
        math::d3 exclusion_low;
        /// exclusion zone for reload: high
        math::d3 exclusion_high;
    };

protected:
    void handlePairList(PairList &pairList) override;

private:
    /// cutoff radius squared
    double m_cutoff2;
    /// is either PMF or updated version
    FunctionPL m_force;
    /// is either PMF or updated version
    FunctionPL m_potential;
    /// exclusion zone for reload: low
    math::d3 m_exclusion_low;
    /// exclusion zone for reload: high
    math::d3 m_exclusion_high;
};



#endif //IBI_FORCEFUNCTOR_H
