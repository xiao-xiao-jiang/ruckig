#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <tuple>
#include <variant>

#include <ruckig/path.hpp>


namespace ruckig {

template<size_t DOFs>
struct Trajectory {
    //! Duration of the synchronized trajectory
    double duration;

    //! Minimum duration of each independent DoF
    std::array<double, DOFs> independent_min_durations;

    //! The path information
    std::variant<ProfileTrajectory<DOFs>, PathTrajectory<DOFs>> _trajectory;

    //! Set of current profiles for each DoF
    std::array<Profile, DOFs> profiles;

    //! Get the output parameter for the given time
    void at_time(double time, std::array<double, DOFs>& new_position, std::array<double, DOFs>& new_velocity, std::array<double, DOFs>& new_acceleration) const {
        std::visit([&](auto&& trajectory) {
            return trajectory.at_time(time, new_position, new_velocity, new_acceleration);
        }, _trajectory);
    }

    //! Get the min/max values of the position for each DoF and the current trajectory
    std::array<PositionExtrema, DOFs> get_position_extrema() const {
        return std::visit([&](auto&& trajectory) {
            return trajectory.get_position_extrema();
        }, _trajectory);
    }
};

} // namespace ruckig
