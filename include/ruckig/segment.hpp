#pragma once

#include <array>


namespace ruckig {

template<size_t DOFs>
struct LinearSegment {
    using Vector = std::array<double, DOFs>;
    
    double length;
    Vector start, end;

    explicit LinearSegment() { }
    explicit LinearSegment(const Vector& start, const Vector& end): start(start), end(end) {
        length = 0.0;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            length += (end[dof] - start[dof]) * (end[dof] - start[dof]);
        }
        length = std::sqrt(length);
    }

    Vector q(double s) const {
        Vector result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = start[dof] + s / length * (end[dof] - start[dof]);
        }
        return result;
    }

    Vector pdq([[maybe_unused]] double s) const {
        Vector result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = (end[dof] - start[dof]) / length;
        }
        return result;
    }

    Vector pddq([[maybe_unused]] double s) const {
        Vector result {};
        return result;
    }

    Vector pdddq([[maybe_unused]] double s) const {
        Vector result {};
        return result;
    }
};


template<size_t DOFs>
class QuarticBlendSegment {
    using Vector = std::array<double, DOFs>;

public:
    double length;
    Vector b, c, e, f;
    Vector lb, lm, rb, rm;

    explicit QuarticBlendSegment() { }
    explicit QuarticBlendSegment(const Vector& lb, const Vector& lm, const Vector& rb, const Vector& rm, double s_mid, double max_diff, double s_abs_max): lb(lb), lm(lm), rb(rb), rm(rm) {
        Vector s_abs;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            s_abs[dof] = std::abs((-16*max_diff)/(3*(lm[dof] - rm[dof])));
        }
         
        const double s_abs_min = std::min<double>(*std::min_element(s_abs.begin(), s_abs.end()), s_abs_max);
        length = 2 * s_abs_min;

        for (size_t dof = 0; dof < DOFs; ++dof) {
            b[dof] = (lm[dof] - rm[dof]) / (16*s_abs_min*s_abs_min*s_abs_min);
            c[dof] = (-lm[dof] + rm[dof]) / (4*s_abs_min*s_abs_min);
            e[dof] = lm[dof];
            f[dof] = lb[dof] + lm[dof]*(s_mid - s_abs_min);
        }
    }

    Vector q(double s) const {
        Vector result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = f[dof] + s * (e[dof] + s * (s * (c[dof] + s * b[dof])));
        }
        return result;
    }

    Vector pdq(double s) const {
        Vector result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = e[dof] + s * (s * (3 * c[dof] + s * 4 * b[dof]));
        }
        return result;
    }

    Vector pddq(double s) const {
        Vector result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = s * (6 * c[dof] + s * 12 * b[dof]);
        }
        return result;
    }

    Vector pdddq(double s) const {
        Vector result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = 6 * c[dof] + s * 24 * b[dof];
        }
        return result;
    }
};

} // namespace ruckig
