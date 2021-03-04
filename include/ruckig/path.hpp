#pragma once

#include <array>
#include <memory>
#include <vector>

#include <ruckig/profile.hpp>


namespace ruckig {

template<size_t DOFs>
struct Segment {
    double s_length;

    virtual std::array<double, DOFs> q(double s) const = 0;
};

} // namespace ruckig
