//
// Created by alex on 8/4/24.
//

#ifndef KOMD_LJ12_6_H
#define KOMD_LJ12_6_H

#include <cstdint>
#include "ForceFunctor.h"
#include "container/SOA.h"

class Site;

class LJ12_6 : public ForceFunctor {
public:
    LJ12_6();
    ~LJ12_6() override = default;

    struct LJ12_6_Force {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA::vec_scatter_t<math::d3> scatter;
        SOA::vec_t<uint64_t> id;
        SOA::vec_t<math::d3> r;
        SOA::vec_t<double> sig;
        SOA::vec_t<double> eps;
        const double cutoff2;
    };

    struct LJ12_6_ForcePair {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA::vec_scatter_t<math::d3> scatter0;
        SOA::vec_scatter_t<math::d3> scatter1;

        SOA::vec_t<uint64_t> id0;
        SOA::vec_t<uint64_t> id1;
        SOA::vec_t<math::d3> r0;
        SOA::vec_t<math::d3> r1;
        SOA::vec_t<double> sig0;
        SOA::vec_t<double> sig1;
        SOA::vec_t<double> eps0;
        SOA::vec_t<double> eps1;
        const double cutoff2;
    };
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    double m_cutoff2;
};


#endif //KOMD_LJ12_6_H
