#pragma once

#include <algorithm>
#include <array>
#include <optional>
#include <variant>
#include <vector>

#include <ruckig/profile.hpp>
#include <ruckig/segment.hpp>


namespace ruckig {

template<size_t DOFs>
struct PathWaypoint {
    using Vector = std::array<double, DOFs>;

    enum class Reference {
        Absolute,
        Relative
    } reference {Reference::Absolute};

    Vector vector;
    std::optional<double> max_blend_distance;

    explicit PathWaypoint(const Vector& vector): vector(vector) { }
    explicit PathWaypoint(const Vector& vector, Reference reference): vector(vector), reference(reference) { }
};


template<size_t DOFs>
class Path {
    using Vector = std::array<double, DOFs>;
    using Segment = std::variant<LinearSegment<DOFs>, QuarticBlendSegment<DOFs>>;

    std::vector<Vector> absolute_waypoints;
    std::vector<LinearSegment<DOFs>> line_segments;

    std::tuple<size_t, double> find_index(double s) const {
        const auto index_ptr = std::upper_bound(cumulative_lengths.begin(), cumulative_lengths.end(), s);
        const size_t index = std::distance(cumulative_lengths.begin(), index_ptr) - 1;
        return {index, s - cumulative_lengths[index]};
    }

public:
    static constexpr size_t degrees_of_freedom {DOFs};
    double length;

    std::vector<Segment> segments;
    std::vector<double> cumulative_lengths;
    
    explicit Path(const Vector& start, const std::vector<PathWaypoint<DOFs>>& waypoints, double max_blend_distance = 0.0) {
        if (waypoints.empty()) {
            // throw std::runtime_error("Path needs at least 2 waypoints as input, but has only " + std::to_string(waypoints.size()) + ".");
        }

        absolute_waypoints.resize(waypoints.size() + 1);
        line_segments.resize(waypoints.size());

        segments.reserve(2*waypoints.size());
        cumulative_lengths.reserve(2*waypoints.size());

        absolute_waypoints[0] = start;
        cumulative_lengths[0] = 0;

        for (size_t i = 0; i < waypoints.size(); ++i) {
            switch (waypoints[i].reference) {
                case PathWaypoint<DOFs>::Reference::Absolute: {
                    absolute_waypoints[i+1] = waypoints[i].vector;
                } break;
                case PathWaypoint<DOFs>::Reference::Relative: {
                    for (size_t dof = 0; dof < DOFs; ++dof) {
                        absolute_waypoints[i+1][dof] = absolute_waypoints[i][dof] + waypoints[i].vector[dof];
                    }
                } break;
            }

            line_segments[i] = LinearSegment(absolute_waypoints[i], absolute_waypoints[i+1]);
        }

        double cumulative_length {0.0};
        cumulative_lengths.emplace_back(cumulative_length);
        for (size_t i = 1; i < line_segments.size(); i += 1) {
            const double blend_distance = waypoints[i].max_blend_distance.value_or(max_blend_distance);
            if (blend_distance > 0.0) {
                auto& left = line_segments[i - 1];
                auto& right = line_segments[i];

                Vector lm, rm;
                for (size_t dof = 0; dof < DOFs; ++dof) {
                    lm[dof] = (left.end[dof] - left.start[dof]) / left.length;
                    rm[dof] = (right.end[dof] - right.start[dof]) / right.length;
                }

                const double s_abs_max = std::min(left.length, right.length) / 2;
                QuarticBlendSegment<DOFs> blend {left.start, lm, right.start, rm, left.length, blend_distance, s_abs_max};
                const double s_abs = blend.length / 2;

                LinearSegment<DOFs> new_left {left.start, left.q(left.length - s_abs)};
                LinearSegment<DOFs> new_right {right.q(s_abs), right.end};

                cumulative_length += new_left.length;
                segments.emplace_back(new_left);
                cumulative_lengths.emplace_back(cumulative_length);

                cumulative_length += blend.length;
                segments.emplace_back(blend);
                cumulative_lengths.emplace_back(cumulative_length);

                right = new_right;

            } else {
                cumulative_length += line_segments[i - 1].length;
                segments.emplace_back(line_segments[i - 1]);
                cumulative_lengths.emplace_back(cumulative_length);
            }
        }

        segments.emplace_back(line_segments.back());
        cumulative_length += line_segments.back().length;
        length = cumulative_length;
    }

    Vector q(double s) const {
        auto [i, s_local] = find_index(s);
        return std::visit([s_local](auto&& segment) { return segment.q(s_local); }, segments[i]);
    }

    Vector pdq(double s) const {
        auto [i, s_local] = find_index(s);
        return std::visit([s_local](auto&& segment) { return segment.pdq(s_local); }, segments[i]);
    }

    Vector pddq(double s) const {
        auto [i, s_local] = find_index(s);
        return std::visit([s_local](auto&& segment) { return segment.pddq(s_local); }, segments[i]);
    }

    Vector pdddq(double s) const {
        auto [i, s_local] = find_index(s);
        return std::visit([s_local](auto&& segment) { return segment.pdddq(s_local); }, segments[i]);
    }

    Vector dq(double s, double ds) const {
        auto [i, s_local] = find_index(s);
        return std::visit([s_local, ds](auto&& segment) {
            Vector result;
            const Vector pdq {segment.pdq(s_local)};
            for (size_t dof = 0; dof < DOFs; ++dof) {
                result[dof] = pdq[dof] * ds;
            }
            return result;
        }, segments[i]);
    }

    Vector ddq(double s, double ds, double dds) const {
        auto [i, s_local] = find_index(s);
        return std::visit([s_local, ds, dds](auto&& segment) {
            Vector result;
            const Vector pdq {segment.pdq(s_local)}, pddq {segment.pddq(s_local)};
            for (size_t dof = 0; dof < DOFs; ++dof) {
                result[dof] = pddq[dof] * ds * ds + pdq[dof] * dds;
            }
            return result;
        }, segments[i]);
    }

    Vector dddq(double s, double ds, double dds, double ddds) const {
        auto [i, s_local] = find_index(s);
        return std::visit([s_local, ds, dds, ddds](auto&& segment) {
            Vector result;
            const Vector pdq {segment.pdq(s_local)}, pddq {segment.pddq(s_local)}, pdddq {segment.pdddq(s_local)};
            for (size_t dof = 0; dof < DOFs; ++dof) {
                result[dof] = 3 * ds * pddq[dof] * dds + ds * ds * ds * pdddq[dof] + pdq[dof] * ddds;
            }
            return result;
        }, segments[i]);
    }
};

} // namespace ruckig
