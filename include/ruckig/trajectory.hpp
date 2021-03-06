#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <tuple>
#include <variant>

#include <ruckig/path.hpp>


namespace ruckig {


enum class CalculationError {
    Working,
};


// For path-based OTG
template<size_t DOFs>
class PathTrajectory {
    using Vector = std::array<double, DOFs>;

    double s0, ds0, dds0;
    double sf, dsf, ddsf;

    Vector p0, v0, a0;
    Vector pf, vf, af;

    std::tuple<double, double, double> time_parametrization(double time) const {
        double s = 0.0;
        double ds = 0.0;
        double dds = 0.0;
        return {s, ds, dds};
    }

public:
    double duration;
    
    Path<DOFs> path;

    explicit PathTrajectory(const Path<DOFs>& path): PathTrajectory(path, {}, {}, {}, {}) { }
    explicit PathTrajectory(const Path<DOFs>& path, Vector v0, Vector a0, Vector vf, Vector af): path(path), s0(0.0), sf(path.length), p0(path.q(s0)), pf(path.q(sf)), v0(v0), vf(vf), a0(a0), af(af) {
        ds0 = 0.0;
        dds0 = 0.0;

        dsf = 0.0;
        ddsf = 0.0;
    }

    CalculationError calculate() {
        return CalculationError::Working;
    }

    void at_time(double time, Vector& new_position, Vector new_velocity, Vector& new_acceleration) const {
        if (time > duration) {
            // Keep constant acceleration from final state
            
        }

        auto [s, ds, dds] = time_parametrization(time);
        new_position = path.q(s);
        new_velocity = path.dq(s, ds);
        new_acceleration = path.ddq(s, ds, dds);
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

    explicit ProfileTrajectory() { }

    CalculationError calculate() {
        return CalculationError::Working;
    }

    void at_time(double time, std::array<double, DOFs>& new_position, std::array<double, DOFs>& new_velocity, std::array<double, DOFs>& new_acceleration) const {
        if (time > duration) {
            // Keep constant acceleration from final state
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
