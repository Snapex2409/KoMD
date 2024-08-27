//
// Created by alex on 7/31/24.
//

#ifndef KOMD_ARRAY_H
#define KOMD_ARRAY_H

#include <type_traits>
#include <cstdint>
#include <array>
#include <cmath>
#include "Kokkos_Core.hpp"

namespace math {
    /**
     * General N-dimensional array template.
     * */
    template<typename T, uint64_t N>
    requires std::is_arithmetic_v<T>
    class Array {
    public:
        /**
         * Constructs Array filled with init value
         * */
        explicit Array(T init) : m_data() {
            for (uint64_t idx = 0; idx < N; idx++) m_data[idx] = init;
        }

        /**
         * Constructs Array filled with 0
         * */
        explicit Array() : m_data() {}

        /// Access to raw pointer
        T* data() { return m_data.data(); }

        Array operator*(const Array& other) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] * other[idx];
            return result;
        }

        Array operator+(const Array& other) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] + other[idx];
            return result;
        }

        Array operator/(const Array& other) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] / other[idx];
            return result;
        }

        Array operator-(const Array& other) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] - other[idx];
            return result;
        }

        Array operator*(T val) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] * val;
            return result;
        }

        Array operator+(T val) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] + val;
            return result;
        }

        Array operator/(T val) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] / val;
            return result;
        }

        Array operator-(T val) const {
            Array<T, N> result;
            for (uint64_t idx = 0; idx < N; idx++) result.m_data[idx] = m_data[idx] - val;
            return result;
        }

        T& operator[](uint64_t idx) {
            return m_data[idx];
        }

        T operator[](uint64_t idx) const {
            return m_data[idx];
        }

        /**
         * Computes dot product with other Array
         * */
        [[maybe_unused]] T dot(const Array& other) const {
            T result = 0;
            for (uint64_t idx = 0; idx < N; idx++) result += m_data[idx] * other.m_data[idx];
            return result;
        }

        /**
         * Computes L1 norm of this array
         * */
        [[maybe_unused]] T L1() const {
            T result = 0;
            for (uint64_t idx = 0; idx < N; idx++) result += std::abs(m_data[idx]);
            return result;
        }

        /**
         * Computes L2 norm of this array
         * */
        [[maybe_unused]] T L2() const {
            T result = 0;
            for (uint64_t idx = 0; idx < N; idx++) result += std::pow(m_data[idx], 2);
            return std::sqrt(result);
        }
    private:
        /// storage of data
        std::array<T, N> m_data;
    };

    /**
     * Specialized 3D array template. Is specialized for Kokkos usage.
     * */
    template<typename T>
    class Array<T, 3>{
    public:
        /**
         * Constructs Array filled with init value
         * */
        KOKKOS_INLINE_FUNCTION explicit Array(T init) : m_data() {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = init;
        }

        /**
         * Constructs Array with provided coordinates
         * */
        KOKKOS_INLINE_FUNCTION constexpr Array(T x, T y, T z) : m_data{x, y, z} {}

        /**
         * Constructs Array filled with 0
         * */
        KOKKOS_INLINE_FUNCTION explicit Array() : m_data() {}

        template<typename O>
        KOKKOS_INLINE_FUNCTION Array& operator=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = other[idx];
            return *this;
        }

        /// Access to raw pointer
        T* data() { return m_data; }

        template<typename O>
        KOKKOS_FUNCTION Array& operator*=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] * other[idx];
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array& operator+=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] + other[idx];
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array& operator/=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] / other[idx];
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array& operator-=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] - other[idx];
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array& operator*=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] * val;
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array& operator+=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] + val;
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array& operator/=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] / val;
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array& operator-=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] - val;
            return *this;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator*(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] * other[idx];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator+(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] + other[idx];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator/(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] / other[idx];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator-(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] - other[idx];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator*(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] * val;
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator+(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] + val;
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator/(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] / val;
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array operator-(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] - val;
            return result;
        }

        KOKKOS_FUNCTION T& operator[](uint64_t idx) {
            return m_data[idx];
        }

        KOKKOS_FUNCTION T operator[](uint64_t idx) const {
            return m_data[idx];
        }

        template<typename O>
        KOKKOS_FUNCTION bool operator==(const Array<O,3>& other) const {
            return other.m_data[0] == m_data[0] && other.m_data[1] == m_data[1] && other.m_data[2] == m_data[2];
        }

        template<typename O>
        KOKKOS_FUNCTION bool operator!=(const Array<O,3>& other) const {
            return other.m_data[0] != m_data[0] || other.m_data[1] != m_data[1] || other.m_data[2] != m_data[2];
        }

        template<typename O>
        KOKKOS_FUNCTION Array<T,3> operator<=(const Array<O,3>& other) const {
            Array<T,3> result {1, 1, 1};
            for (int dim = 0; dim < 3; dim++) result[dim] = m_data[dim] <= other.m_data[dim];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array<T,3> operator<(const Array<O,3>& other) const {
            Array<T,3> result {1, 1, 1};
            for (int dim = 0; dim < 3; dim++) result[dim] = m_data[dim] < other.m_data[dim];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array<T,3> operator>=(const Array<O,3>& other) const {
            Array<T,3> result {1, 1, 1};
            for (int dim = 0; dim < 3; dim++) result[dim] = m_data[dim] >= other.m_data[dim];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION Array<T,3> operator>(const Array<O,3>& other) const {
            Array<T,3> result {1, 1, 1};
            for (int dim = 0; dim < 3; dim++) result[dim] = m_data[dim] > other.m_data[dim];
            return result;
        }

        template<typename O>
        KOKKOS_FUNCTION bool operator==(O val) const {
            return val == m_data[0] && val == m_data[1] && val == m_data[2];
        }

        template<typename O>
        KOKKOS_FUNCTION bool operator!=(O val) const {
            return val != m_data[0] || val != m_data[1] || val != m_data[2];
        }

        /**
         * Computes dot product of this and another array
         * */
        template<typename O>
        [[maybe_unused]] KOKKOS_FUNCTION T dot(const Array<O,3>& other) const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += m_data[idx] * other.m_data[idx];
            return result;
        }

        /**
         * Computes L1 norm of this array
         * */
        [[maybe_unused]] KOKKOS_FUNCTION T L1() const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += Kokkos::abs(m_data[idx]);
            return result;
        }

        /**
         * Computes L2 norm of this array
         * */
        [[maybe_unused]] KOKKOS_FUNCTION T L2() const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += Kokkos::pow(m_data[idx], 2);
            return Kokkos::sqrt(result);
        }

        /**
         * Computes sum of elements
         * */
        KOKKOS_FUNCTION T sum() const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += m_data[idx];
            return result;
        }

        /**
         * Computes product of elements
         * */
        KOKKOS_FUNCTION T product() const {
            T result = 1;
            for (uint64_t idx = 0; idx < 3; idx++) result *= m_data[idx];
            return result;
        }

        /// @returns value at dim 0
        KOKKOS_INLINE_FUNCTION T& x() { return m_data[0]; }
        /// @returns value at dim 1
        KOKKOS_INLINE_FUNCTION T& y() { return m_data[1]; }
        /// @returns value at dim 2
        KOKKOS_INLINE_FUNCTION T& z() { return m_data[2]; }
        /// @returns value at dim 0
        KOKKOS_INLINE_FUNCTION T x() const { return m_data[0]; }
        /// @returns value at dim 1
        KOKKOS_INLINE_FUNCTION T y() const { return m_data[1]; }
        /// @returns value at dim 2
        KOKKOS_INLINE_FUNCTION T z() const { return m_data[2]; }

    private:
        /// storage of data
        T m_data[3];
    };

    using i3 = Array<int32_t, 3>;
    using ui3 = Array<uint32_t, 3>;
    using l3 = Array<int64_t, 3>;
    using ul3 = Array<uint64_t, 3>;
    using d3 = Array<double, 3>;
    using f3 = Array<float, 3>;

    [[maybe_unused]] KOKKOS_INLINE_FUNCTION i3 ceil(const f3 &vec) {
        i3 result {0, 0, 0};
        result[0] = Kokkos::ceil(vec[0]);
        result[1] = Kokkos::ceil(vec[1]);
        result[2] = Kokkos::ceil(vec[2]);
        return result;
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION i3 floor(const f3& vec) {
        i3 result {0, 0, 0};
        result[0] = Kokkos::floor(vec[0]);
        result[1] = Kokkos::floor(vec[1]);
        result[2] = Kokkos::floor(vec[2]);
        return result;
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION l3 ceil(const d3& vec) {
        l3 result {0, 0, 0};
        result[0] = Kokkos::ceil(vec[0]);
        result[1] = Kokkos::ceil(vec[1]);
        result[2] = Kokkos::ceil(vec[2]);
        return result;
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION l3 floor(const d3& vec) {
        l3 result {0, 0, 0};
        result[0] = Kokkos::floor(vec[0]);
        result[1] = Kokkos::floor(vec[1]);
        result[2] = Kokkos::floor(vec[2]);
        return result;
    }

    [[maybe_unused]] KOKKOS_INLINE_FUNCTION ui3 uceil(const f3 &vec) {
        ui3 result {0, 0, 0};
        result[0] = Kokkos::ceil(vec[0]);
        result[1] = Kokkos::ceil(vec[1]);
        result[2] = Kokkos::ceil(vec[2]);
        return result;
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION ui3 ufloor(const f3& vec) {
        ui3 result {0, 0, 0};
        result[0] = Kokkos::floor(vec[0]);
        result[1] = Kokkos::floor(vec[1]);
        result[2] = Kokkos::floor(vec[2]);
        return result;
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION ul3 uceil(const d3& vec) {
        ul3 result {0, 0, 0};
        result[0] = Kokkos::ceil(vec[0]);
        result[1] = Kokkos::ceil(vec[1]);
        result[2] = Kokkos::ceil(vec[2]);
        return result;
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION ul3 ufloor(const d3& vec) {
        ul3 result {0, 0, 0};
        result[0] = Kokkos::floor(vec[0]);
        result[1] = Kokkos::floor(vec[1]);
        result[2] = Kokkos::floor(vec[2]);
        return result;
    }

    [[maybe_unused]] KOKKOS_INLINE_FUNCTION d3 max(const d3 &a, const d3 &b) {
        return {Kokkos::max(a.x(), b.x()), Kokkos::max(a.y(), b.y()), Kokkos::max(a.z(), b.z())};
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION d3 min(const d3 &a, const d3 &b) {
        return {Kokkos::min(a.x(), b.x()), Kokkos::min(a.y(), b.y()), Kokkos::min(a.z(), b.z())};
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION d3 max(const d3 &a, double v) {
        return {Kokkos::max(a.x(), v), Kokkos::max(a.y(), v), Kokkos::max(a.z(), v)};
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION d3 min(const d3 &a, double v) {
        return {Kokkos::min(a.x(), v), Kokkos::min(a.y(), v), Kokkos::min(a.z(), v)};
    }
    [[maybe_unused]] KOKKOS_INLINE_FUNCTION d3 abs(const d3 &vec) {
        return {Kokkos::abs(vec.x()), Kokkos::abs(vec.y()), Kokkos::abs(vec.z())};
    }

} // math

#endif //KOMD_ARRAY_H
