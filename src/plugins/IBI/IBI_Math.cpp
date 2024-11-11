//
// Created by alex on 9/20/24.
//
#include "IBI_Math.h"
#include <fstream>
#include <vector>
#include <IO/Logging.h>

FunctionPL::FunctionPL(uint64_t buffer_size, double def_low, double def_high): default_value_lower(def_low), default_value_upper(def_high) {
    x_values = KW::vec_t<double>("Function PL x", buffer_size);
    y_values = KW::vec_t<double>("Function PL y", buffer_size);
}

void FunctionPL::setXValues(const KW::vec_t<double>& x) {
    for (int idx = 0; idx < x_values.size(); ++idx) x_values[idx] = x[idx];
}

void FunctionPL::setYValues(const KW::vec_t<double>& y) {
    for (int idx = 0; idx < y_values.size(); ++idx) y_values[idx] = y[idx];
}

void FunctionPL::derivative(FunctionPL& target) const {
    target.setXValues(x_values);

    const auto size = y_values.size();
    // forward finite diff at start
    target.y_values[0] = (y_values[1] - y_values[0]) / (x_values[1] - x_values[0]);
    // backward finite diff at end
    target.y_values[size - 1] = (y_values[size - 1] - y_values[size - 2]) / (x_values[size - 1] - x_values[size - 2]);
    // central finite diff for rest
    for (int idx = 1; idx < size - 1; idx++) {
        target.y_values[idx] = (y_values[idx + 1] - y_values[idx - 1]) / ((x_values[idx + 1] - x_values[idx - 1]));
    }
}

void FunctionPL::read(const std::string &path) {
    std::ifstream file{path};
    if (!file) {
        Log::io->error() << "Could not read the file" << std::endl;
        throw std::runtime_error("Failed to open file.");
    }

    std::vector<double> x_tmp {};
    std::vector<double> y_tmp {};
    double n1, n2;
    while (file >> n1 >> n2) {
        x_tmp.push_back(n1);
        y_tmp.push_back(n2);
    }
    file.close();

    if (x_values.size() != x_tmp.size()) Kokkos::resize(x_values, x_tmp.size());
    if (y_values.size() != y_tmp.size()) Kokkos::resize(y_values, y_tmp.size());

    for (int idx = 0; idx < x_tmp.size(); idx++) {
        x_values[idx] = x_tmp[idx];
        y_values[idx] = y_tmp[idx];
    }
}

void FunctionPL::write(const std::string &path) {
    std::ofstream file{path};
    if (!file) {
        Log::io->error() << "Could not read the file" << std::endl;
        throw std::runtime_error("Failed to open file.");
    }

    bool begin = true;
    for (int idx = 0; idx < x_values.size(); idx++) {
        if (begin) begin = false;
        else file << "\n";

        file << x_values[idx] << "\t" << y_values[idx];
    }
    file.close();
}
