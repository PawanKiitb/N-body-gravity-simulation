# Documentation for CPU implementation of N-body simulation

## Overview
So the implementation of the N-body simulation on the CPU is a two step process.
- Firstly, we **compute the forces** acting on each body due to the gravitational attraction of all other bodies and store them in the acceleration vector of each body.
- Secondly, we **update the positions and velocities** of each body based on the computed forces.

## File structure
The implementation is organized into the following files:
- **common.hpp:** Contains all the important constants, data structures, and utility functions used across the implementation.
- **particle.hpp:** Contains the definition of a three-dimensional vector class, and the definition of the Particle class with its properties and methods.
- **force_calculation.hpp:** Contains the implementations of several methods for calculating the forces between particles.
- **time_integration_methods.cpp:** Contains the implementations of several methods for updating the positions and velocities of particles based on the computed forces.
- **n-body.cpp:** Is the main file that orchestrates the simulation by calling the appropriate functions for force calculation and time integration.

## File Descriptions

### common.hpp
Included important libraries:
- `<iostream>`: For input and output operations.
-  `<cmath>`: For mathematical functions.
- `<string>`: For storing names of particles.
- `<vector>`: Useful container
- `<fstream>`: For storing data to files.

Important constants:
- `PI` = 3.14159265358979323846: Value of pi
- `G` = 39.4784176044: Gravitational constant
- `AU` = 1.495978707e11: Astronomical unit in meters
- `SOLAR_MASS`   = 1.98847e30: Mass of the sun in kilograms
- `SOLAR_YEAR`   = 31557600: Number of seconds in a year

### particle.hpp
This file contains the definition of a three-dimensional vector class, and the definition of the Particle class with its properties and methods.

**Vector3D class:**
- Represents a three-dimensional vector with x, y, and z components.
- Provides methods for vector addition, subtraction, scalar multiplication, and computing the magnitude of the vector
- Provides methods for computing the dot product and cross product of two vectors.

**Particle class:**
- Represents a particle with mass, position, velocity, acceleration, and a name(prefrably the name of the celestial body or a unique identifier).

### force_calculation.hpp
