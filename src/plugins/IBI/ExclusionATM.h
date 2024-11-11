//
// Created by alex on 11/11/24.
//

#ifndef EXCLUSIONATM_H
#define EXCLUSIONATM_H

#include "potentials/ForceFunctor3B.h"

class ExclusionATM : public ForceFunctor3B {
public:
    ExclusionATM();
    ~ExclusionATM() override = default;

    /// Kernel for device
    struct ExATM_Force {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of TripleList::triplets
        KW::nvec_t<int, 3> triplets;
        /// shallow copy of TripleList::triple_offsets
        KW::nvec_t<math::d3, 3> triple_offsets;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// Energy parameter
        const double nu;
        /// cutoff radius squared
        const double cutoff2;
        /// exclusion zone for reload: low
        math::d3 exclusion_low;
        /// exclusion zone for reload: high
        math::d3 exclusion_high;

        static KOKKOS_INLINE_FUNCTION double comp_force(const double nu, const double r_ab, const double r_ac, const double r_bc) {
            const double r_ab2 = r_ab  * r_ab;
            const double r_ab3 = r_ab2 * r_ab;
            const double r_ab4 = r_ab2 * r_ab2;
            //const double r_ab5 = r_ab3 * r_ab2;
            const double r_ab6 = r_ab3 * r_ab3;

            const double r_ac2 = r_ac  * r_ac;
            const double r_ac3 = r_ac2 * r_ac;
            //const double r_ac4 = r_ac2 * r_ac2;
            const double r_ac5 = r_ac3 * r_ac2;
            //const double r_ac6 = r_ac3 * r_ac3;

            const double r_bc2 = r_bc  * r_bc;
            const double r_bc3 = r_bc2 * r_bc;
            //const double r_bc4 = r_bc2 * r_bc2;
            const double r_bc5 = r_bc3 * r_bc2;
            //const double r_bc6 = r_bc3 * r_bc3;

            return ((3. * nu) / (8. * r_ab)) *
                   (
                       - (8.)/(r_ab4 * r_ac3 * r_bc3)
                       - (1.)/(r_ac5 * r_bc5)
                       + (5. * r_ac)/(r_ab6 * r_bc5)
                       + (5. * r_bc)/(r_ab6 * r_ac5)
                       - (1.)/(r_ab2 * r_ac3 * r_bc5)
                       - (1.)/(r_ab2 * r_ac5 * r_bc3)
                       - (3.)/(r_ab4 * r_ac * r_bc5)
                       - (3.)/(r_ab4 * r_ac5 * r_bc)
                       - (5.)/(r_ab6 * r_ac * r_bc3)
                       - (5.)/(r_ab6 * r_ac3 * r_bc)
                       + (6.)/(r_ab4 * r_ac3 * r_bc3)
                   );
        }
    };

protected:
    void handleTripleList(TripleList &tripleList) override;

private:
    /// Energy parameter
    double m_nu;
    /// cutoff radius squared
    double m_cutoff2;
    /// exclusion zone for reload: low
    math::d3 m_exclusion_low;
    /// exclusion zone for reload: high
    math::d3 m_exclusion_high;
};



#endif //EXCLUSIONATM_H
