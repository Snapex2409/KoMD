//
// Created by alex on 7/31/24.
//

#ifndef KOMD_SITE_H
#define KOMD_SITE_H

#include "math/Array.h"
#include "util/defaults.h"

/**
 * Singular Particle using LJ potential
 * */
class Site {
public:
    explicit Site(double epsilon = Defaults::epsilon, double sigma = Defaults::sigma, double mass = Defaults::mass, const math::d3& r = Defaults::r, const math::d3& v = Defaults::v);

    /// @returns epsilon
    [[maybe_unused, nodiscard]] double getEpsilon() const;
    /// @returns sigma
    [[maybe_unused, nodiscard]] double getSigma() const;
    /// @returns mass
    [[maybe_unused, nodiscard]] double getMass() const;
    /// @returns position
    [[maybe_unused, nodiscard]] math::d3& r_arr();
    /// @returns force
    [[maybe_unused, nodiscard]] math::d3& f_arr();
    /// @returns velocity
    [[maybe_unused, nodiscard]] math::d3& v_arr();
    /// @returns position
    [[maybe_unused, nodiscard]] const math::d3& r_arr() const;
    /// @returns force
    [[maybe_unused, nodiscard]] const math::d3& f_arr() const;
    /// @returns velocity
    [[maybe_unused, nodiscard]] const math::d3& v_arr() const;
private:
    /// epsilon value
    double m_epsilon;
    /// sigma value
    double m_sigma;
    /// mass
    double m_m;
    /// position
    math::d3 m_r;
    /// force
    math::d3 m_f;
    /// velocity
    math::d3 m_v;
};


#endif //KOMD_SITE_H
