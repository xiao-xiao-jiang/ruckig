#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <tuple>
#include <variant>

#include <ruckig/path.hpp>
#include <ruckig/input_parameter.hpp>
#include <ruckig/position.hpp>
#include <ruckig/velocity.hpp>
#include <ruckig/block.hpp>
#include <ruckig/segment.hpp>


namespace ruckig {

// For path-based OTG
template<size_t DOFs>
class PathTrajectory {
    using Vector = std::array<double, DOFs>;
    constexpr static double eps {1e-10};

    double s0, ds0, dds0;
    double sf, dsf, ddsf;

    Vector p0, v0, a0;
    Vector pf, vf, af;

    Profile main_profile;

    std::vector<double> cumulative_times;

    std::tuple<double, double, double> time_parametrization(double time) const {
        const auto [new_p, new_v, new_a] = main_profile.state_at_time(time);

        const auto segment = std::get<LinearSegment<DOFs>>(path.segments[0]);
        const double scale = (main_profile.pf - main_profile.p[0]) / segment.length;
        double s = new_p / scale;
        double ds = new_v / scale;
        double dds = new_a / scale;
        return {s, ds, dds};
    }

public:
    double duration;
    std::array<double, DOFs> independent_min_durations;
    
    Path<DOFs> path;

    explicit PathTrajectory(const Path<DOFs>& path, Vector v0, Vector a0, Vector vf, Vector af): path(path), s0(0.0), sf(path.length), p0(path.q(0.0)), pf(path.q(path.length)), v0(v0), vf(vf), a0(a0), af(af) {}

    bool validate_boundary() {
        auto pdq_s0 = path.pdq(s0);
        auto pddq_s0 = path.pddq(s0);
        auto pdq_sf = path.pdq(sf);
        auto pddq_sf = path.pddq(sf);

        // Calculate from first DoF
        ds0 = v0[0] / pdq_s0[0];
        dds0 = (a0[0] - pddq_s0[0] * ds0 * ds0) / pdq_s0[0];

        dsf = vf[0] / pdq_sf[0];
        ddsf = (af[0] - pddq_sf[0] * dsf * dsf) / pdq_sf[0];

        for (size_t dof = 1; dof < DOFs; ++dof) {
            const double ds0_dof = v0[dof] / pdq_s0[dof];
            const double dds0_dof = (a0[dof] - pddq_s0[dof] * ds0 * ds0) / pdq_s0[dof];

            const double dsf_dof = vf[dof] / pdq_sf[dof];
            const double ddsf_dof = (af[dof] - pddq_sf[dof] * dsf * dsf) / pdq_sf[dof];

            if (std::abs(ds0 - ds0_dof) > eps || std::abs(dds0 - dds0_dof) > eps || std::abs(dsf - dsf_dof) > eps || std::abs(ddsf - ddsf_dof) > eps) {
                return false;
            }
        }

        return true;
    }

    template<bool throw_error>
    CalculationResult calculate(const ruckig::InputParameter<DOFs>& input, double delta_time) {
        if (path.segments.size() != 1 || !std::holds_alternative<LinearSegment<DOFs>>(path.segments[0])) {
            std::cout << "WIP: Currently not supported." << std::endl;
        }

        std::array<Profile, DOFs> profiles;
        std::array<Block, DOFs> blocks;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            Block block;
            Position1 step1 {p0[dof], v0[dof], a0[dof], pf[dof], vf[dof], af[dof], input.max_velocity[dof], -input.max_velocity[dof], input.max_acceleration[dof], -input.max_acceleration[dof], input.max_jerk[dof]};
            step1.get_profile(profiles[dof], blocks[dof]);
        }

        double t_sync;
        int limiting_dof;
        Block::synchronize(blocks, std::nullopt, duration, limiting_dof, profiles, false, 0.001); 

        if (limiting_dof < 0) {
            // TODO: Use step2
        }
        
        main_profile = profiles[limiting_dof];
        duration = main_profile.t_sum[6];
        return CalculationResult::Working;
    }

    void at_time(double time, Vector& new_position, Vector new_velocity, Vector& new_acceleration) const {
        if (time > duration) {
            // Keep constant acceleration from final state
            for (size_t dof = 0; dof < DOFs; ++dof) {
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(time - duration, pf[dof], vf[dof], af[dof], 0);
            }
            return;
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
    constexpr static double eps {std::numeric_limits<double>::epsilon()};

    double duration;
    std::array<double, DOFs> independent_min_durations;

    //! Set of current profiles for each DoF
    std::array<Profile, DOFs> profiles;

    explicit ProfileTrajectory() { }

    template<bool throw_error, bool return_error_at_maximal_duration>
    CalculationResult calculate(const InputParameter<DOFs>& inp, double delta_time) {
        std::array<Block, DOFs> blocks;
        std::array<double, DOFs> p0s, v0s, a0s; // Starting point of profiles without brake trajectory
        std::array<double, DOFs> inp_min_velocity, inp_min_acceleration;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            auto& p = profiles[dof];

            if (!inp.enabled[dof]) {
                p.pf = inp.current_position[dof];
                p.vf = inp.current_velocity[dof];
                p.af = inp.current_acceleration[dof];
                p.t_sum[6] = 0.0;
                continue;
            }

            inp_min_velocity[dof] = inp.min_velocity ? inp.min_velocity.value()[dof] : -inp.max_velocity[dof];
            inp_min_acceleration[dof] = inp.min_acceleration ? inp.min_acceleration.value()[dof] : -inp.max_acceleration[dof];

            // Calculate brake (if input exceeds or will exceed limits)
            switch (inp.interface) {
                case Interface::Position: {
                    Brake::get_position_brake_trajectory(inp.current_velocity[dof], inp.current_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof], p.t_brakes, p.j_brakes);
                } break;
                case Interface::Velocity: {
                    Brake::get_velocity_brake_trajectory(inp.current_acceleration[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof], p.t_brakes, p.j_brakes);
                } break;
            }

            p.t_brake = p.t_brakes[0] + p.t_brakes[1];
            p0s[dof] = inp.current_position[dof];
            v0s[dof] = inp.current_velocity[dof];
            a0s[dof] = inp.current_acceleration[dof];

            // Integrate brake pre-trajectory
            for (size_t i = 0; p.t_brakes[i] > 0 && i < 2; ++i) {
                p.p_brakes[i] = p0s[dof];
                p.v_brakes[i] = v0s[dof];
                p.a_brakes[i] = a0s[dof];
                std::tie(p0s[dof], v0s[dof], a0s[dof]) = Profile::integrate(p.t_brakes[i], p0s[dof], v0s[dof], a0s[dof], p.j_brakes[i]);
            }

            bool found_profile;
            switch (inp.interface) {
                case Interface::Position: {
                    Position1 step1 {p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_profile = step1.get_profile(p, blocks[dof]);
                } break;
                case Interface::Velocity: {
                    Velocity1 step1 {p0s[dof], v0s[dof], a0s[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_profile = step1.get_profile(p, blocks[dof]);
                } break;
            }

            if (!found_profile) {
                if constexpr (throw_error) {
                    throw std::runtime_error("[ruckig] error in step 1, dof: " + std::to_string(dof) + " input: " + inp.to_string());
                }
                return CalculationResult::ErrorExecutionTimeCalculation;
            }

            independent_min_durations[dof] = blocks[dof].t_min;
        }

        int limiting_dof; // The DoF that doesn't need step 2
        const bool discrete_duration = (inp.duration_discretization == DurationDiscretization::Discrete);
        const bool found_synchronization = Block::synchronize(blocks, inp.minimum_duration, duration, limiting_dof, profiles, discrete_duration, delta_time);
        if (!found_synchronization) {
            if constexpr (throw_error) {
                throw std::runtime_error("[ruckig] error in time synchronization: " + std::to_string(duration));
            }
            return CalculationResult::ErrorSynchronizationCalculation;
        }

        if constexpr (return_error_at_maximal_duration) {
            if (duration > 7.6e3) {
                return CalculationResult::ErrorTrajectoryDuration;
            }
        }

        if (duration > 0.0 && inp.synchronization != Synchronization::None) {
            for (size_t dof = 0; dof < DOFs; ++dof) {
                if (!inp.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                Profile& p = profiles[dof];
                const double t_profile = duration - p.t_brake.value_or(0.0);

                if (inp.synchronization == Synchronization::TimeIfNecessary && std::abs(inp.target_velocity[dof]) < eps && std::abs(inp.target_acceleration[dof]) < eps) {
                    p = blocks[dof].p_min;
                    continue;
                }

                // Check if the final time corresponds to an extremal profile calculated in step 1
                if (std::abs(t_profile - blocks[dof].t_min) < eps) {
                    p = blocks[dof].p_min;
                    continue;
                } else if (blocks[dof].a && std::abs(t_profile - blocks[dof].a->right) < eps) {
                    p = blocks[dof].a->profile;
                    continue;
                } else if (blocks[dof].b && std::abs(t_profile - blocks[dof].b->right) < eps) {
                    p = blocks[dof].b->profile;
                    continue;
                }

                bool found_time_synchronization;
                switch (inp.interface) {
                    case Interface::Position: {
                        Position2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                        found_time_synchronization = step2.get_profile(p);
                    } break;
                    case Interface::Velocity: {
                        Velocity2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                        found_time_synchronization = step2.get_profile(p);
                    } break;
                }
                if (!found_time_synchronization) {
                    if constexpr (throw_error) {
                        throw std::runtime_error("[ruckig] error in step 2 in dof: " + std::to_string(dof) + " for t sync: " + std::to_string(duration) + " input: " + inp.to_string());
                    }
                    return CalculationResult::ErrorSynchronizationCalculation;
                }
                // std::cout << dof << " profile step2: " << p.to_string() << std::endl;
            }

        } else if (inp.synchronization == Synchronization::None) {
            for (size_t dof = 0; dof < DOFs; ++dof) {
                if (!inp.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                profiles[dof] = blocks[dof].p_min;
            }
        }

        return CalculationResult::Working;
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

            std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = p.state_at_time(t_diff);
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

    template<class T>
    void set_data(const T& trajectory) {
        data = trajectory;
        duration = trajectory.duration;
        independent_min_durations = trajectory.independent_min_durations;
    }

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
