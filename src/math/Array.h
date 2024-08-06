//
// Created by alex on 7/31/24.
//

#ifndef KOMD_ARRAY_H
#define KOMD_ARRAY_H

#include <type_traits>
#include <cstdint>
#include <array>
#include <cmath>

namespace math {
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

        [[maybe_unused]] T dot(const Array& other) const {
            T result = 0;
            for (uint64_t idx = 0; idx < N; idx++) result += m_data[idx] * other.m_data[idx];
            return result;
        }

        [[maybe_unused]] T L1() const {
            T result = 0;
            for (uint64_t idx = 0; idx < N; idx++) result += std::abs(m_data[idx]);
            return result;
        }

        [[maybe_unused]] T L2() const {
            T result = 0;
            for (uint64_t idx = 0; idx < N; idx++) result += std::pow(m_data[idx], 2);
            return std::sqrt(result);
        }
    private:
        /// storage of data
        std::array<T, N> m_data;
    };

    template<typename T>
    class Array<T, 3>{
    public:
        /**
         * Constructs Array filled with init value
         * */
        explicit Array(T init) : m_data() {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = init;
        }

        /**
         * Constructs Array with provided coordinates
         * */
        constexpr Array(T x, T y, T z) : m_data{x, y, z} {}

        /**
         * Constructs Array filled with 0
         * */
        explicit Array() : m_data() {}

        /// Access to raw pointer
        T* data() { return m_data.data(); }

        template<typename O>
        Array& operator*=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] * other[idx];
            return *this;
        }

        template<typename O>
        Array& operator+=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] + other[idx];
            return *this;
        }

        template<typename O>
        Array& operator/=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] / other[idx];
            return *this;
        }

        template<typename O>
        Array& operator-=(const Array<O,3>& other) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] - other[idx];
            return *this;
        }

        template<typename O>
        Array& operator*=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] * val;
            return *this;
        }

        template<typename O>
        Array& operator+=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] + val;
            return *this;
        }

        template<typename O>
        Array& operator/=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] / val;
            return *this;
        }

        template<typename O>
        Array& operator-=(O val) {
            for (uint64_t idx = 0; idx < 3; idx++) m_data[idx] = m_data[idx] - val;
            return *this;
        }

        template<typename O>
        Array operator*(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] * other[idx];
            return result;
        }

        template<typename O>
        Array operator+(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] + other[idx];
            return result;
        }

        template<typename O>
        Array operator/(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] / other[idx];
            return result;
        }

        template<typename O>
        Array operator-(const Array<O,3>& other) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] - other[idx];
            return result;
        }

        template<typename O>
        Array operator*(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] * val;
            return result;
        }

        template<typename O>
        Array operator+(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] + val;
            return result;
        }

        template<typename O>
        Array operator/(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] / val;
            return result;
        }

        template<typename O>
        Array operator-(O val) const {
            Array<T, 3> result;
            for (uint64_t idx = 0; idx < 3; idx++) result.m_data[idx] = m_data[idx] - val;
            return result;
        }

        T& operator[](uint64_t idx) {
            return m_data[idx];
        }

        T operator[](uint64_t idx) const {
            return m_data[idx];
        }

        template<typename O>
        bool operator==(const Array<O,3>& other) {
            return other.m_data[0] == m_data[0] && other.m_data[1] == m_data[1] && other.m_data[2] == m_data[2];
        }

        template<typename O>
        bool operator!=(const Array<O,3>& other) {
            return other.m_data[0] != m_data[0] || other.m_data[1] != m_data[1] || other.m_data[2] != m_data[2];
        }

        template<typename O>
        [[maybe_unused]] T dot(const Array<O,3>& other) const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += m_data[idx] * other.m_data[idx];
            return result;
        }

        [[maybe_unused]] T L1() const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += std::abs(m_data[idx]);
            return result;
        }

        [[maybe_unused]] T L2() const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += std::pow(m_data[idx], 2);
            return std::sqrt(result);
        }

        T sum() const {
            T result = 0;
            for (uint64_t idx = 0; idx < 3; idx++) result += m_data[idx];
            return result;
        }

        T product() const {
            T result = 1;
            for (uint64_t idx = 0; idx < 3; idx++) result *= m_data[idx];
            return result;
        }

        inline T& x() { return m_data[0]; }
        inline T& y() { return m_data[1]; }
        inline T& z() { return m_data[2]; }
        inline T x() const { return m_data[0]; }
        inline T y() const { return m_data[1]; }
        inline T z() const { return m_data[2]; }

    private:
        /// storage of data
        std::array<T, 3> m_data;
    };

    using i3 = Array<int32_t, 3>;
    using ui3 = Array<uint32_t, 3>;
    using l3 = Array<int64_t, 3>;
    using ul3 = Array<uint64_t, 3>;
    using d3 = Array<double, 3>;
    using f3 = Array<float, 3>;

    [[maybe_unused]] i3 ceil(const f3 &vec);
    [[maybe_unused]] i3 floor(const f3& vec);
    [[maybe_unused]] l3 ceil(const d3& vec);
    [[maybe_unused]] l3 floor(const d3& vec);

    [[maybe_unused]] ui3 uceil(const f3 &vec);
    [[maybe_unused]] ui3 ufloor(const f3& vec);
    [[maybe_unused]] ul3 uceil(const d3& vec);
    [[maybe_unused]] ul3 ufloor(const d3& vec);

    [[maybe_unused]] d3 max(const d3 &a, const d3 &b);
    [[maybe_unused]] d3 min(const d3 &a, const d3 &b);
    [[maybe_unused]] d3 max(const d3 &a, double v);
    [[maybe_unused]] d3 min(const d3 &a, double v);
    [[maybe_unused]] d3 abs(const d3 &vec);

} // math

#endif //KOMD_ARRAY_H
