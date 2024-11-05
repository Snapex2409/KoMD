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
    KOKKOS_FUNCTION static double evaluateAt(double x, double def_low, double def_high, const KW::vec_t<double>& x_values, const KW::vec_t<double>& y_values);
    void derivative(FunctionPL& target) const;

    template<typename T>
    FunctionPL& operator*=(T val) { for (int idx = 0; idx < y_values.size(); idx++) y_values[idx] *= val; return *this; }
    template<typename T>
    FunctionPL& operator+=(T val) { for (int idx = 0; idx < y_values.size(); idx++) y_values[idx] *= val; return *this; }

    /// read from file
    void read(const std::string& path);
    /// write to file
    void write(const std::string& path);
private:
    /// Interpolates between two points
    static double linearInterpolation(int a, int b, double x, const KW::vec_t<double>& x_values, const KW::vec_t<double>& y_values);
    /// Returns the index of the next large x node
    static int findUpperNode(double x, const KW::vec_t<double>& x_values);
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
