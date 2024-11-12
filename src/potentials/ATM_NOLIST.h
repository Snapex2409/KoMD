//
// Created by alex on 11/12/24.
//

#ifndef ATM_NOLIST_H
#define ATM_NOLIST_H

#include "ForceFunctor3B.h"
#include "util/Kokkos_Wrapper.h"

class ATM_NOLIST : public ForceFunctor3B {
public:
    ATM_NOLIST();
    ~ATM_NOLIST() override = default;

    /// "Kernel" for device, now is called from launcher kernel down below
    struct ATM_NOLIST_Force {
        KOKKOS_FUNCTION void operator()(uint64_t s_idx_0, uint64_t s_idx_1, uint64_t s_idx_2, const math::d3& offset0, const math::d3& offset1, const math::d3& offset2) const;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// Energy parameter
        double nu;
        /// cutoff radius squared
        double cutoff2;
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

    struct ATM_NOLIST_CellTripleData {
        math::d3 shift0;
        math::d3 shift1;
        math::d3 shift2;
        KW::vec_t<uint64_t> indices0;
        KW::vec_t<uint64_t> indices1;
        KW::vec_t<uint64_t> indices2;
    };

    struct ATM_NOLIST_CellData {
        KW::vec_t<uint64_t> indices;
    };

    struct ATM_NOLISTTriple_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const {
            // this is a global index across multiple cell triplets -> need to find correct cell triple first
            int triple_idx = -1;
            for (int i = 0; i < stored_cell_triplets; i++) {
                if (idx < triple_counts[i] + triple_counts_accumulated[i]) {
                    triple_idx = i;
                    break;
                }
            }
            if (triple_idx == -1) return;

            // we found our triple -> need to get x, y, z coords
            const int local_idx = idx - triple_counts_accumulated[triple_idx];
            const int i = local_idx % triple_dims[triple_idx][0];
            const int j = (local_idx / triple_dims[triple_idx][0]) % triple_dims[triple_idx][1]; // needs to same dim, but this time we just divide...
            const int k = local_idx / (triple_dims[triple_idx][0] * triple_dims[triple_idx][1]);

            // execute force
            force_kernel(cell_triplets[triple_idx].indices0(i),
                         cell_triplets[triple_idx].indices1(j),
                         cell_triplets[triple_idx].indices2(k),
                         cell_triplets[triple_idx].shift0,
                         cell_triplets[triple_idx].shift1,
                         cell_triplets[triple_idx].shift2);
        }

        static constexpr int MAX_TRIPLETS = 256;
        KW::Array<ATM_NOLIST_CellTripleData, MAX_TRIPLETS> cell_triplets;
        KW::Array<int, MAX_TRIPLETS> triple_counts;
        KW::Array<int, MAX_TRIPLETS> triple_counts_accumulated;
        KW::Array<KW::Array<int, 3>, MAX_TRIPLETS> triple_dims;
        int stored_cell_triplets = 0;
        ATM_NOLIST_Force force_kernel {};
    };

    struct ATM_NOLISTSingle_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const {
            // this is a global index across multiple cell triplets -> need to find correct cell triple first
            int triple_idx = -1;
            for (int i = 0; i < stored_cells; i++) {
                if (idx < counts[i] + counts_accumulated[i]) {
                    triple_idx = i;
                    break;
                }
            }
            if (triple_idx == -1) return;

            // we found our pair -> need to get x, y and z coords
            const int local_idx = idx - counts_accumulated[triple_idx];
            int k_idx = -1;
            for (uint64_t pos = 0; pos < tri_sums.size(); pos++) {
                if (tri_sums[pos] > local_idx) {
                    k_idx = pos;
                    break;
                }
            }

            const int k = k_idx+2;
            int tmp_idx = local_idx;
            if (k_idx > 0) tmp_idx -= tri_sums[k_idx-1];
            int j = Kokkos::floor((-1 + Kokkos::sqrt(1 + 8 * tmp_idx)) / 2.0);
            const int i = tmp_idx - j * (j+1) / 2;
            j += 1;

            // execute force
            force_kernel(cells[triple_idx].indices(i),
                         cells[triple_idx].indices(j),
                         cells[triple_idx].indices(k),
                         {0, 0, 0},
                         {0, 0, 0},
                         {0, 0, 0});
        }

        static constexpr int MAX_TRIPLETS = 64;
        KW::Array<ATM_NOLIST_CellData, MAX_TRIPLETS> cells;
        KW::Array<int, MAX_TRIPLETS> counts;
        KW::Array<int, MAX_TRIPLETS> counts_accumulated;
        KW::vec_t<int> tri_sums; // idx 0 has S_1, this are the sums of the first k triangle numbers
        int stored_cells = 0;
        ATM_NOLIST_Force force_kernel {};
    };

protected:
    void handleTripleList(TripleList &tripleList) override;

private:
    /// Energy parameter
    double m_nu;
    /// cutoff radius squared
    double m_cutoff2;
    /// domain low
    math::d3 m_low;
    /// domain high
    math::d3 m_high;
};



#endif //ATM_NOLIST_H
