//
// Created by alex on 7/31/24.
//

#include "Site.h"

Site::Site(double epsilon, double sigma, double mass, const math::d3 &r, const math::d3 &v) : m_epsilon(epsilon), m_sigma(sigma), m_m(mass), m_r(r), m_v(v) { }

[[maybe_unused, nodiscard]] double Site::getEpsilon() const {
    return m_epsilon;
}
[[maybe_unused, nodiscard]] double Site::getSigma() const {
    return m_sigma;
}
[[maybe_unused, nodiscard]] double Site::getMass() const {
    return m_m;
}
[[maybe_unused, nodiscard]] math::d3 &Site::r_arr() {
    return m_r;
}
[[maybe_unused, nodiscard]] math::d3 &Site::f_arr() {
    return m_f;
}
[[maybe_unused, nodiscard]] math::d3 &Site::v_arr() {
    return m_v;
}
[[maybe_unused, nodiscard]] const math::d3 &Site::r_arr() const {
    return m_r;
}
[[maybe_unused, nodiscard]] const math::d3 &Site::f_arr() const {
    return m_f;
}
[[maybe_unused, nodiscard]] const math::d3 &Site::v_arr() const {
    return m_v;
}
