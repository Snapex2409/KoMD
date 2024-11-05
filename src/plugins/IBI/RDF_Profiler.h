//
// Created by alex on 11/4/24.
//

#ifndef RDF_PROFILER_H
#define RDF_PROFILER_H
#include "util/Kokkos_Wrapper.h"

/**
* Specialized RDF profiler for IBI purposes. Exposes more API.
* */
class RDF_Profiler {
public:
    RDF_Profiler();
    /// creates a new sample for RDF
    void measure();
    /// sets buffers back to 0
    void reset();
    /// writes current normalized rdf into provided buffer (must have correct size)
    void getRDF(KW::vec_t<double>& buffer) const;
    /// returns rdf buffer size
    [[nodiscard]] uint64_t getBufferSize() const;
    /// writes r nodes into buffer (must have correct size)
    void getRNodes(KW::vec_t<double>& buffer) const;
    /// returns cref to nodes
    [[nodiscard]] const KW::vec_t<double>& getRNodes() const { return m_nodes; }

    /// Kernel for device
    struct IBI_RDF_Kernel {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        /// shallow copy of CoM positions
        KW::vec_t<math::d3> positions;
        /// shallow copy of rdf bins
        KW::vec_t<double> bins;
        /// maximum measure distance
        const double max_r;
        /// bin size
        const double delta_r;
        /// number of total bins
        const uint64_t num_bins;
        /// domain size
        const math::d3 domain_size;
    };
private:
    /**
     * @param r distance r in U(r)
     * @param delta_r size of one bin
     * @param bins number of bins
     * */
    static KOKKOS_INLINE_FUNCTION uint64_t getBin(double r, double delta_r, uint64_t bins) {
        if (r < 0) return bins-1;

        auto bin = static_cast<uint64_t>(r / delta_r);
        return Kokkos::clamp(bin, 0UL, bins - 1);
    }

    /// maximum measure distance
    double m_max_r;
    /// bin size
    double m_delta_r;
    /// rdf buffer
    KW::vec_t<double> m_bins;
    /// number of samples
    double m_samples;
    /// com buffer
    KW::vec_t<math::d3> m_com;
    /// measure node positions
    KW::vec_t<double> m_nodes;
};

#endif //RDF_PROFILER_H
