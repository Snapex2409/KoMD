//
// Created by alex on 10/15/24.
//

#include "IBI.h"
#include "Registry.h"
#include "IO/Logging.h"

IBI::IBI() : Plugin("IBI"),
reference_rdf(Registry::instance->configuration()->IBI_bins, 0.0, 1.0),
reference_potential(Registry::instance->configuration()->IBI_bins, 1e+6, 0.0),
update_function(Registry::instance->configuration()->IBI_bins, 0, 0)
{
    auto config = Registry::instance->configuration();
    alpha = config->IBI_alpha;
    T = config->temperature;
    steps_equilibration = config->IBI_steps_equil;
    steps_measurement = config->IBI_steps_measure;

    update_function.setXValues(profiler.getRNodes());

    convergenceCheck.init(
        config->IBI_conv_threshold,
        config->IBI_conv_mode,
        config->IBI_conv_stop,
        config->IBI_conv_window);
}

void IBI::post_forces() {
    current_steps++;

        switch (ibi_phase) {
            case FIRST_INIT: {
                if (current_steps < steps_measurement) { profiler.measure(); break; }

                ibi_phase = EQUILIBRATE;
                current_steps = 0;
                ibi_iteration = 0;

                // create rdf_0
                profiler.getRDF(reference_rdf.getYValues());
                profiler.getRNodes(reference_rdf.getXValues());
                reference_rdf.write(createFilepath("rdf"));

                profiler.reset();

                // set up remaining fields
                initializePotentialValues();

                // replace force function to IBI implementation
                force_functor = std::make_shared<IBI_ForceFunctor>();
                Registry::instance->forceFunctors().clear();
                Registry::instance->forceFunctors3b().clear();
                Registry::instance->forceFunctors().push_back(force_functor);

                force_functor->getPotentialFunction().setXValues(reference_potential.getXValues());
                force_functor->getPotentialFunction().setYValues(reference_potential.getYValues());
                force_functor->getPotentialFunction().setLowerDefault(reference_potential.getLowerDefault());
                force_functor->getPotentialFunction().setUpperDefault(reference_potential.getUpperDefault());

                derivativeOfPotential();

                Log::simulation->info() << "[IBI] Transitioning from first initialization to equilibration with PMF" << std::endl;
                break;
            }
            case EQUILIBRATE: {
                if (current_steps < steps_equilibration) break;
                ibi_phase = MEASURE;
                current_steps = 0;
                // nothing else to do in this phase
                Log::simulation->info() << "[IBI] Transitioning from equilibration to measurement" << std::endl;
                break;
            }
            case MEASURE: {
                if (current_steps < steps_measurement) { profiler.measure(); break; }

                ibi_phase = EQUILIBRATE;
                current_steps = 0;

                // update pot -> derivative
                addPotentialCorrection();
                ibi_iteration++;
                force_functor->getPotentialFunction().write(createFilepath("pot"));
                derivativeOfPotential();
                writeRDF();

                auto conv = convergenceCheck(reference_rdf, profiler, update_function.getYValues());
                Log::simulation->info() << "[IBI] Convergence: target_reached=" << conv.first << " value=" << conv.second << std::endl;
                const bool should_stop = convergenceCheck.shouldStop();
                if (should_stop) {
                    Log::simulation->info() << "[IBI] Convergence: stopping criterion reached. Stopping now." << std::endl;
                    Registry::instance->simulation()->stop();
                    break;
                }

                Log::simulation->info() << "[IBI] Transitioning from measurement to equilibration" << std::endl;
                profiler.reset();
                break;
            }
        }
}

void IBI::finalize() {
    force_functor->getPotentialFunction().write("final_pot.txt");

    std::stringstream ss;
    convergenceCheck.logValues(ss);
    Log::simulation->info() << "[IBI] Convergence steps: " << ss.str() << std::endl;
}

void IBI::initializePotentialValues() {
    reference_potential.setXValues(reference_rdf.getXValues());
    reference_potential.setYValues(reference_rdf.getYValues());

    auto& u0 = reference_potential.getYValues();
    for (int idx = 0; idx < u0.size(); idx++) u0[idx] = -T * Kokkos::log(u0[idx]);

    // extrapolate function values (linearly) for which log was undefined in range of [0, r_min)
    // first find r_min
    int r_min_idx = -1;
    for (int idx = 0; idx < reference_potential.getXValues().size(); idx++) {
        if (Kokkos::isfinite(reference_potential.getYValues()[idx])) { r_min_idx = idx; break; }
    }
    // there are inf values -> need to do something
    // r_min_idx has the first position with a valid y value
    if (r_min_idx != -1) {
        const double dy = reference_potential.getYValues()[r_min_idx + 1] - reference_potential.getYValues()[r_min_idx];
        double y_val = reference_potential.getYValues()[r_min_idx];
        for (int idx = r_min_idx; idx >= 0; idx--) {
            reference_potential.getYValues()[idx] = y_val;
            y_val -= dy;
        }
    }

    reference_potential.write("pot_0.txt");
}

void IBI::addPotentialCorrection() {
    profiler.getRDF(update_function.getYValues());
    auto& pot = force_functor->getPotentialFunction().getYValues();

    for (int idx = 0; idx < pot.size(); idx++) {
        const double update = alpha * T * Kokkos::log(update_function.getYValues()[idx] / reference_rdf.getYValues()[idx]);
        if (Kokkos::isfinite(update)) {
            pot[idx] -= update;
            update_function.getYValues()[idx] = update;
        }
        else {
            update_function.getYValues()[idx] = 0;
        }
    }
    update_function.write(createFilepath("update"));

    // extrapolate function values (linearly) for which log was undefined in range of [0, r_min)
    // first find r_min
    int r_min_idx = -1;
    for (int idx = 0; idx < reference_rdf.getXValues().size(); idx++) {
        if (reference_rdf.getYValues()[idx] != 0.0) { r_min_idx = idx; break; }
    }
    // r_min_idx has the first position with a valid y value
    if (r_min_idx != -1) {
        const double dy = pot[r_min_idx + 1] - pot[r_min_idx];
        double y_val = pot[r_min_idx];
        for (int idx = r_min_idx; idx >= 0; idx--) {
            pot[idx] = y_val;
            y_val -= dy;
        }
    }
}

void IBI::derivativeOfPotential() {
    force_functor->getPotentialFunction().derivative(force_functor->getForceFunction());
    force_functor->getForceFunction().setLowerDefault(1e+6);
    force_functor->getForceFunction().setUpperDefault(0.0);
    force_functor->getForceFunction() *= -1;
    force_functor->getForceFunction().write(createFilepath("force"));
}

void IBI::writeRDF() {
    profiler.getRDF(update_function.getYValues());
    update_function.write(createFilepath("rdf"));
}

std::string IBI::createFilepath(const std::string &prefix) const {
    std::stringstream path;
    path << prefix << "_" << ibi_iteration << ".txt";
    return path.str();
}

void IBI::Convergence::init(double th, const std::string &mode_str, const std::string &stop_str, int window) {
    if (mode_str == "l2") mode = L2;
    else mode = INTEGRAL;
    threshold = th;
    if (stop_str == "worse") stopping_mode = ON_WORSE;
    else stopping_mode = WINDOW;
    window_size = window;
}

std::pair<bool, double> IBI::Convergence::integral(const FunctionPL &ref, RDF_Profiler &profiler, KW::vec_t<double>& g_i) {
    const auto& x = ref.getXValues();
    const auto& g_0 = ref.getYValues();
    profiler.getRDF(g_i);

    // integrate diff
    double diff = 0;
    for (int idx = 0; idx < x.size()-1; idx++) {
        const double dx = x[idx+1] - x[idx];
        const double lower = Kokkos::abs(g_i[idx] - g_0[idx]);
        const double upper = Kokkos::abs(g_i[idx+1] - g_0[idx+1]);
        diff += dx * (upper + lower) / 2.0;
    }

    // integrate sum
    double sum = 0;
    for (int idx = 0; idx < x.size()-1; idx++) {
        const double dx = x[idx+1] - x[idx];
        const double lower = std::abs(g_i[idx] + g_0[idx]);
        const double upper = std::abs(g_i[idx+1] + g_0[idx+1]);
        sum += dx * (upper + lower) / 2.0;
    }

    if (sum == 0) sum = 1e-15;
    const double conv_value = 1.0 - (diff / sum);
    conv_values.push_back(conv_value);
    return {conv_value <= threshold, conv_value};
}

std::pair<bool, double> IBI::Convergence::l2(const FunctionPL &ref, RDF_Profiler &profiler, KW::vec_t<double>& g_i) {
    const auto& g_0 = ref.getYValues();
    profiler.getRDF(g_i);

    double vec_norm = 0.0;
    for (int i = 0; i < g_0.size(); ++i) {
        vec_norm += Kokkos::pow(g_0[i] - g_i[i], 2.0);
    }

    conv_values.push_back(vec_norm);
    return {vec_norm <= threshold, vec_norm};
}

std::pair<bool, double> IBI::Convergence::operator()(const FunctionPL &ref, RDF_Profiler &profiler, KW::vec_t<double>& g_i) {
    if (mode == L2) return l2(ref, profiler, g_i);
    else if (mode == INTEGRAL) return integral(ref, profiler, g_i);
    else throw std::runtime_error("Unknown convergence method");
}

void IBI::Convergence::logValues(std::ostream &ostream) {
    for (const double conv : conv_values) {
        ostream << conv << " ";
    }
}

bool IBI::Convergence::shouldStop() const {
    const auto steps = conv_values.size();

    if (stopping_mode == ON_WORSE) {
        if (steps < 2) return false;

        if (conv_values[steps-2] > threshold) return false;
        if (conv_values[steps-1] <= threshold) return false;
        return true;
    }
    else if (stopping_mode == WINDOW) {
        if (steps < window_size) return false;
        bool all_leq = true;
        for (int i = 0; i < window_size; i++) {
            all_leq &= conv_values[steps - 1 - i] <= threshold;
        }

        return all_leq;
    }
    throw std::runtime_error("Unknown stopping method");
}
