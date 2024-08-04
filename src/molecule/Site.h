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

    [[maybe_unused, nodiscard]] double getEpsilon() const;
    [[maybe_unused, nodiscard]] double getSigma() const;
    [[maybe_unused, nodiscard]] double getMass() const;
    [[maybe_unused, nodiscard]] math::d3& r_arr();
    [[maybe_unused, nodiscard]] math::d3& f_arr();
    [[maybe_unused, nodiscard]] math::d3& fold_arr();
    [[maybe_unused, nodiscard]] math::d3& v_arr();
    [[maybe_unused, nodiscard]] const math::d3& r_arr() const;
    [[maybe_unused, nodiscard]] const math::d3& f_arr() const;
    [[maybe_unused, nodiscard]] const math::d3& fold_arr() const;
    [[maybe_unused, nodiscard]] const math::d3& v_arr() const;
private:
    double m_epsilon;
    double m_sigma;
    double m_m;
    math::d3 m_r;
    math::d3 m_f;
    math::d3 m_fold;
    math::d3 m_v;
};


#endif //KOMD_SITE_H
