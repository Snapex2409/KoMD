//
// Created by alex on 11/12/24.
//

#include "ATM_NOLIST.h"

#include "Registry.h"
#include "container/Cell.h"
#include "container/LinkedCells.h"
#include "util/constants.h"

ATM_NOLIST::ATM_NOLIST() :
    m_nu(Registry::instance->configuration()->energy_3b),
    m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)),
    m_low(Registry::instance->configuration()->domainLow),
    m_high(Registry::instance->configuration()->domainHigh) {
    m_nu = m_nu * 1e-9 * Constants::conv_J_Ei; // convert external unit to internal
}

void ATM_NOLIST::handleTripleList(TripleList &) {
    // when this class is used, then the passed triple list is useless (its empty)
    auto container = Registry::instance->moleculeContainer();
    // create single force kernel
    ATM_NOLIST_Force force_kernel {p_soa.f(), p_soa.r(), m_nu, m_cutoff2};

    // create cell triple kernels
    std::vector<ATM_NOLISTTriple_Kernel> triple_kernels;
    ATM_NOLISTTriple_Kernel triple_kernel;
    triple_kernel.force_kernel = force_kernel;
    triple_kernel.stored_cell_triplets = 0;

    for (auto it = container->iteratorTripletC08(); it->isValid(); ++(*it)) {
        Cell &cell0 = it->cell0();
        Cell &cell1 = it->cell1();
        Cell &cell2 = it->cell2();
        const math::d3 shift0 = it->getCell0Shift();
        const math::d3 shift1 = it->getCell1Shift();
        const math::d3 shift2 = it->getCell2Shift();

        if (triple_kernel.stored_cell_triplets == ATM_NOLISTTriple_Kernel::MAX_TRIPLETS) {
            triple_kernels.push_back(triple_kernel);
            triple_kernel = ATM_NOLISTTriple_Kernel();
            triple_kernel.force_kernel = force_kernel;
            triple_kernel.stored_cell_triplets = 0;
        }

        const int write_idx = triple_kernel.stored_cell_triplets;
        // cell data
        triple_kernel.cell_triplets[write_idx].indices0 = cell0.indices();
        triple_kernel.cell_triplets[write_idx].indices1 = cell1.indices();
        triple_kernel.cell_triplets[write_idx].indices2 = cell2.indices();
        triple_kernel.cell_triplets[write_idx].shift0 = shift0;
        triple_kernel.cell_triplets[write_idx].shift1 = shift1;
        triple_kernel.cell_triplets[write_idx].shift2 = shift2;

        // meta data
        triple_kernel.triple_counts[write_idx] = cell0.getNumIndices() * cell1.getNumIndices() * cell2.getNumIndices();
        int count_acc = 0;
        if (write_idx > 0) count_acc = triple_kernel.triple_counts_accumulated[write_idx - 1] + triple_kernel.
                                       triple_counts[write_idx - 1];
        triple_kernel.triple_counts_accumulated[write_idx] = count_acc;
        triple_kernel.triple_dims[write_idx][0] = cell0.getNumIndices();
        triple_kernel.triple_dims[write_idx][1] = cell1.getNumIndices();
        triple_kernel.triple_dims[write_idx][2] = cell2.getNumIndices();
        triple_kernel.stored_cell_triplets++;
    }
    if (triple_kernel.stored_cell_triplets != 0) triple_kernels.push_back(triple_kernel);

    // create single cell kernels
    std::vector<ATM_NOLISTSingle_Kernel> single_kernels;
    ATM_NOLISTSingle_Kernel single_kernel;
    single_kernel.force_kernel = force_kernel;
    single_kernel.stored_cells = 0;

    unsigned long single_cell_max_n = 0;
    for (auto it = container->iteratorCell(); it->isValid(); ++(*it)) {
        Cell &cell = it->cell();
        const auto check_low = cell.low() < m_low;
        const auto check_high = cell.high() > m_high;
        const auto is_halo = check_low - check_high; // will be 1 if is halo low, 0 if domain, -1 if is halo high
        if (is_halo != 0) continue;

        if (single_kernel.stored_cells == ATM_NOLISTSingle_Kernel::MAX_TRIPLETS) {
            single_kernels.push_back(single_kernel);
            single_kernel = ATM_NOLISTSingle_Kernel();
            single_kernel.force_kernel = force_kernel;
            single_kernel.stored_cells = 0;
        }

        const int write_idx = single_kernel.stored_cells;
        // cell data
        single_kernel.cells[write_idx].indices = cell.indices();

        // meta data
        const auto n = cell.getNumIndices();
        single_cell_max_n = std::max(single_cell_max_n, n);
        single_kernel.counts[write_idx] = n * (n - 1) * (n - 2) / 6;
        int count_acc = 0;
        if (write_idx > 0) count_acc = single_kernel.counts_accumulated[write_idx - 1] + single_kernel.counts[
                                           write_idx - 1];
        single_kernel.counts_accumulated[write_idx] = count_acc;
        single_kernel.stored_cells++;
    }
    if (single_kernel.stored_cells != 0) single_kernels.push_back(single_kernel);

    KW::vec_t<int> tri_sums = KW::vec_t<int>("ATM tri sums", single_cell_max_n);
    for (int i = 1; i < single_cell_max_n + 1; i++) {
        tri_sums[i - 1] = i * (i + 1) * (i + 2) / 6;
    }
    for (auto &kernel: single_kernels) {
        kernel.tri_sums = tri_sums;
    }

    // run kernels
    for (auto &kernel: triple_kernels) Kokkos::parallel_for("ATM - Triple",
                                                            kernel.triple_counts_accumulated[
                                                                kernel.stored_cell_triplets - 1] + kernel.triple_counts[
                                                                kernel.stored_cell_triplets - 1], kernel);
    for (auto &kernel: single_kernels) Kokkos::parallel_for("ATM - Single",
                                                            kernel.counts_accumulated[kernel.stored_cells - 1] + kernel.
                                                            counts[kernel.stored_cells - 1], kernel);
}

void ATM_NOLIST::ATM_NOLIST_Force::operator()(uint64_t s_idx_0, uint64_t s_idx_1, uint64_t s_idx_2,
                                              const math::d3 &shift0,
                                              const math::d3 &shift1,
                                              const math::d3 &shift2) const {
    const math::d3 r0 = r[s_idx_0] + shift0;
    const math::d3 r1 = r[s_idx_1] + shift1;
    const math::d3 r2 = r[s_idx_2] + shift2;

    const math::d3 dr0_1 = r0 - r1;
    const math::d3 dr0_2 = r0 - r2;
    const math::d3 dr1_2 = r1 - r2;

    const double dr0_1_2 = dr0_1.dot(dr0_1);
    const double dr0_2_2 = dr0_2.dot(dr0_2);
    const double dr1_2_2 = dr1_2.dot(dr1_2);
    // symmetrical condition
    if (dr0_1_2 > cutoff2 || dr0_2_2 > cutoff2 || dr1_2_2 > cutoff2) return;

    // we are within range
    // compute forces
    const double r0_1 = Kokkos::sqrt(dr0_1_2);
    const double r0_2 = Kokkos::sqrt(dr0_2_2);
    const double r1_2 = Kokkos::sqrt(dr1_2_2);

    const double dVdR01 = comp_force(nu, r0_1, r0_2, r1_2);
    const double dVdR02 = comp_force(nu, r0_2, r0_1, r1_2);

    const math::d3 f0 = dr0_1 * dVdR01 + dr0_2 * dVdR02;
    const math::d3 f1 = dr0_1 * (-dVdR01) + dr0_2 * dVdR02;
    const math::d3 f2 = dr0_1 * (-dVdR01) + dr0_2 * (-dVdR02);

    if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], f0);
    if (shift1 == 0) Kokkos::atomic_add(&f[s_idx_1], f1);
    if (shift2 == 0) Kokkos::atomic_add(&f[s_idx_2], f2);
}
