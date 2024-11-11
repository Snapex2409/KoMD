//
// Created by alex on 10/15/24.
//

#ifndef IBI_H
#define IBI_H

#include "IBI_ForceFunctor.h"
#include "IBI_Math.h"
#include "RDF_Profiler.h"
#include "plugins/Plugin.h"
#include "math/Array.h"
#include "util/Kokkos_Wrapper.h"

class IBI : public Plugin {
public:
    IBI();

    ~IBI() override = default;

    void post_forces() override;

    void finalize() override;

private:
    /// Implements U(r)_0 = -T^*ln(g(r)_0)
    void initializePotentialValues();
    /// Implements U(r)_{i+1} =U(r)_i - alpha*T^*ln(g(r)_i/g(r)_*)
    void addPotentialCorrection();
    /// Implements F(r) = - d/dr U(r)
    void derivativeOfPotential();
    /// Writes the total average RDF
    void writeRDF();
    /// Generates file paths
    [[nodiscard]] std::string createFilepath(const std::string& prefix) const;

    /// g*(r)
    FunctionPL reference_rdf;
    /// V0(r)
    FunctionPL reference_potential;
    /// measuring tool
    RDF_Profiler profiler;
    /// pair handler to comp forces
    std::shared_ptr<IBI_ForceFunctor<IBI_Default>> force_functor;
    /// alpha for step size
    double alpha = 0.0;
    /// Target simulation temperature
    double T = 0.0;
    /// IBI iteration counter
    int ibi_iteration = 0;
    /// number of equilibration steps
    int steps_equilibration = 1e+6;
    /// number of rdf measurement steps
    int steps_measurement = 1e+5;
    /// counter for elapsed simulation steps of current IBI phase
    int current_steps = 0;
    /// single run: current ibi phase
    enum { FIRST_INIT, EQUILIBRATE, MEASURE } ibi_phase = FIRST_INIT;
    /// buffer to store tmp rdf
    FunctionPL update_function;

    /**
     * Provides methods for convergence checking
     * */
    struct Convergence {
        void init(double threshold, const std::string& mode_str = "l2", const std::string& stop_str = "worse", int window = 10);

        /**
         * Computes integrals over r_i and r_0... see https://arxiv.org/pdf/1410.1853
         * */
        std::pair<bool, double> integral(const FunctionPL& ref, RDF_Profiler& profiler, KW::vec_t<double>& rdf_buffer);

        /**
         * Computes simple l2 norm over functions
         * */
        std::pair<bool, double> l2(const FunctionPL& ref, RDF_Profiler& profiler, KW::vec_t<double>& rdf_buffer);

        std::pair<bool, double> operator()(const FunctionPL &ref, RDF_Profiler& profiler, KW::vec_t<double>& rdf_buffer);

        void logValues(std::ostream& ostream);

        bool shouldStop() const;

    private:
        /// Convergence values of each iteration, yes this is intended (should be CPU only)
        std::vector<double> conv_values;
        /// mode selection
        enum {INTEGRAL, L2} mode = L2;
        /// conv threshold
        double threshold;
        /// early stopping method
        enum {ON_WORSE, WINDOW} stopping_mode = ON_WORSE;
        /// window size for window stopping method
        int window_size;

    } convergenceCheck;
};



#endif //IBI_H
