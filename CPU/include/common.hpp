/***** Important Includes *****/


#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

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
}

