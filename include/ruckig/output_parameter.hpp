#pragma once

#include <array>
#include <iomanip>
#include <optional>
#include <sstream>

#include <ruckig/input_parameter.hpp>
#include <ruckig/path.hpp>
#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Output type of the OTG
template<size_t DOFs>
struct OutputParameter {
    static constexpr size_t degrees_of_freedom {DOFs};
    
    std::array<double, DOFs> new_position, new_velocity, new_acceleration;

    //! Was a new trajectory calculation performed in the last cycle?
    bool new_calculation {false};
    double calculation_duration; // [µs]

    //! Current trajectory
    Trajectory<DOFs> trajectory;

    //! Current time on trajectory
    double time;

    //! Current type
    Type type;
};

} // namespace ruckig
