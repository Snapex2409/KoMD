//
// Created by alex on 9/20/24.
//

#ifndef IBI_MATH_H
#define IBI_MATH_H

#include <string>
#include "util/Kokkos_Wrapper.h"

/**
 * Piecewise linear function
 * */
class FunctionPL {
public:
    explicit FunctionPL(uint64_t buffer_size, double def_low = 0.0, double def_high = 0.0);
    /// x must have same size as local buffers
    void setXValues(const KW::vec_t<double>& x);
    /// y must have same size as local buffers
    void setYValues(const KW::vec_t<double>& y);

    [[nodiscard]] const KW::vec_t<double>& getXValues() const { return x_values; };
    KW::vec_t<double>& getXValues() { return x_values; };
    [[nodiscard]] const KW::vec_t<double>& getYValues() const { return y_values; };
    KW::vec_t<double>& getYValues() { return y_values; };
    [[nodiscard]] double getLowerDefault() const { return default_value_lower; }
    [[nodiscard]] double getUpperDefault() const { return default_value_upper; }
    void setLowerDefault(double value) { default_value_lower = value; }
    void setUpperDefault(double value) { default_value_upper = value; }
    KOKKOS_INLINE_FUNCTION static double evaluateAt(double x, double def_low, double def_high, const KW::vec_t<double>& x_values, const KW::vec_t<double>& y_values) {
        if (x > x_values[x_values.size() - 1]) return def_high;
        if (x < x_values[0]) return def_low;

        //find between which 2 values
        const auto idx_upper = findUpperNode(x, x_values);
        const auto idx_lower = idx_upper -1;

        //interpolate
        const auto fx = linearInterpolation(idx_lower, idx_upper, x, x_values, y_values);
        return fx;
    }
    void derivative(FunctionPL& target) const;

    template<typename T>
    FunctionPL& operator*=(T val) { for (int idx = 0; idx < y_values.size(); idx++) y_values[idx] *= val; return *this; }
    template<typename T>
    FunctionPL& operator+=(T val) { for (int idx = 0; idx < y_values.size(); idx++) y_values[idx] *= val; return *this; }

    /// read from file
    void read(const std::string& path);
    /// write to file
    void write(const std::string& path);

    /// Interpolates between two points
    KOKKOS_INLINE_FUNCTION static double linearInterpolation(int a, int b, double x, const KW::vec_t<double>& x_values, const KW::vec_t<double>& y_values) {
        const auto ya = y_values[a];
        const auto yb = y_values[b];
        const auto xa = x_values[a];
        const auto xb = x_values[b];

        return ya + (yb - ya) / (xb - xa) * (x - xa);
    }
    /// Returns the index of the next large x node
    KOKKOS_INLINE_FUNCTION static int findUpperNode(double x, const KW::vec_t<double>& x_values) {
        int i = 0;
        const int i_max = x_values.size();
        while(x > x_values[i]) {
            i++;
            if(i == i_max) break;
        }
        return i;
    }
private:
    /// function value nodes
    KW::vec_t<double> x_values;
    /// function values
    KW::vec_t<double> y_values;
    /// default value lower
    double default_value_lower;
    /// default value upper
    double default_value_upper;
};
#endif //IBI_MATH_H
