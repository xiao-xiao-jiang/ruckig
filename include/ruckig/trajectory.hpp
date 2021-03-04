#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <tuple>
#include <variant>

#include <ruckig/path.hpp>


namespace ruckig {

// For path-based OTG
template<size_t DOFs>
struct PathTrajectory {
    double duration;

    std::vector<double> cumulative_lengths {1};
    std::vector<std::shared_ptr<Segment<DOFs>>> segments {1};

    double s_length() {
        return cumulative_lengths.back();
    }

    void at_time(double time, std::array<double, DOFs>& new_position, std::array<double, DOFs>& new_velocity, std::array<double, DOFs>& new_acceleration) const {
    }

    std::array<PositionExtrema, DOFs> get_position_extrema() const {
        std::array<PositionExtrema, DOFs> result;
        return result;
    }
};


// For waypoint-based OTG
template<size_t DOFs>
struct ProfileTrajectory {
    double duration;

    //! Set of current profiles for each DoF
    std::array<Profile, DOFs> profiles;

    void at_time(double time, std::array<double, DOFs>& new_position, std::array<double, DOFs>& new_velocity, std::array<double, DOFs>& new_acceleration) const {
        if (time > duration) {
            // Keep constant acceleration
            for (size_t dof = 0; dof < DOFs; ++dof) {
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(time - duration, profiles[dof].pf, profiles[dof].vf, profiles[dof].af, 0);
            }
            return;
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            const Profile& p = profiles[dof];

            double t_diff = time;
            if (p.t_brake) {
                if (t_diff < p.t_brake.value()) {
                    const size_t index = (t_diff < p.t_brakes[0]) ? 0 : 1;
                    if (index > 0) {
                        t_diff -= p.t_brakes[index - 1];
                    }

                    std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff, p.p_brakes[index], p.v_brakes[index], p.a_brakes[index], p.j_brakes[index]);
                    continue;
                } else {
                    t_diff -= p.t_brake.value();
                }
            }

            // Non-time synchronization
            if (t_diff >= p.t_sum[6]) {
                // Keep constant acceleration
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff - p.t_sum[6], p.pf, p.vf, p.af, 0);
                continue;
            }

            const auto index_ptr = std::upper_bound(p.t_sum.begin(), p.t_sum.end(), t_diff);
            const size_t index = std::distance(p.t_sum.begin(), index_ptr);

            if (index > 0) {
                t_diff -= p.t_sum[index - 1];
            }

            std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff, p.p[index], p.v[index], p.a[index], p.j[index]);
        }
    }

    std::array<PositionExtrema, DOFs> get_position_extrema() const {
        std::array<PositionExtrema, DOFs> result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = profiles[dof].get_position_extrema();
        }
        return result;
    }
};


template<size_t DOFs>
struct Trajectory {
    //! Duration of the synchronized trajectory
    double duration;

    //! Minimum duration of each independent DoF
    std::array<double, DOFs> independent_min_durations;

    //! The path information
    std::variant<ProfileTrajectory<DOFs>, PathTrajectory<DOFs>> data;

    //! Get the output parameter for the given time
    void at_time(double time, std::array<double, DOFs>& new_position, std::array<double, DOFs>& new_velocity, std::array<double, DOFs>& new_acceleration) const {
        std::visit([&](auto&& trajectory) {
            return trajectory.at_time(time, new_position, new_velocity, new_acceleration);
        }, data);
    }

    //! Get the min/max values of the position for each DoF and the current trajectory
    std::array<PositionExtrema, DOFs> get_position_extrema() const {
        return std::visit([&](auto&& trajectory) {
            return trajectory.get_position_extrema();
        }, data);
    }
};

} // namespace ruckig
