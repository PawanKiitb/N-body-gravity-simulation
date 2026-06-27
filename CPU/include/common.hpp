/***** Important Includes *****/

#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>

// Alias for clock
using Clock = std::chrono::high_resolution_clock;

/***** Important Constants *****/
/*
We aren't going to use SI units
We'll be using solar systems units

Unit System:
Distance : Astronomical Units (AU)
Mass     : Solar Masses (M_sun)
Time     : Years (yr)

In this system:
G = 4π²
*/
namespace constant {
    constexpr double PI = 3.14159265358979323846; // Value of Pi
    constexpr double G = 4.0 * PI * PI;
    constexpr double AU = 1.495978707e11;      // meters
    constexpr double SOLAR_MASS = 1.98847e30;  // kg
    constexpr double SOLAR_YEAR = 365.25 * 24 * 3600; // seconds

    // Other used constants
    constexpr size_t MAX_NODES = 1000000;
    constexpr size_t MAX_PARTICLES = 1000000;
    constexpr int leaf_capacity = 8; // we'll find the optimal value for this later
}

/****** Timing Method *******/
struct TimingStatistics {
    std::string name = "";
    double total_time = 0.0;
    size_t calls = 0;
    double average_time() const {
        return calls > 0 ? total_time / calls : 0.0;
    }
    void reset() {
        total_time = 0.0;
        calls = 0;
    }
};