#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cub/cub.cuh>


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
    constexpr int MAX_NODES = 1000000;
    constexpr int MAX_PARTICLES = 1000000;
    constexpr int LEAF_CAPACITY = 8;
    constexpr int MAX_DEPTH = 10;
    constexpr double EPSILON = 1e-9;

    // CUDA specific constants
    constexpr int BLOCK_SIZE = 256; // Number of threads per block
    constexpr int TILE_SIZE = BLOCK_SIZE; // Tile size for shared memory optimization
}