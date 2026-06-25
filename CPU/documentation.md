# Documentation for CPU Implementation of N-body Simulation

## Overview

This project implements an object-oriented gravitational N-body simulator running on the CPU. The simulation follows Newton's law of universal gravitation and updates the motion of particles using different numerical integration schemes.

The simulation consists of two independent stages performed every timestep:

1. **Force Calculation**

   * Compute the gravitational acceleration acting on every particle due to all other particles.
   * Different algorithms can be used for this stage (Naive (O(N^2)), Barnes-Hut, etc.).

2. **Time Integration**

   * Update particle positions and velocities using the computed accelerations.
   * Different numerical integrators (Euler, Euler-Cromer, Leapfrog, Velocity Verlet, RK4, etc.) can be selected independently of the force calculation method.

This separation allows any force calculation algorithm to be paired with any numerical integrator.

---

# Project Structure

```
common.hpp

particle.hpp

force_calculator.hpp
force_calculator_naive.hpp
force_calculator_naive.cpp

integrator.hpp
integrator_euler.hpp
integrator_euler.cpp
integrator_euler_cromer.hpp
integrator_euler_cromer.cpp
integrator_leapfrog.hpp
integrator_leapfrog.cpp
integrator_velocity_verlet.hpp
integrator_velocity_verlet.cpp
integrator_rk4.hpp
integrator_rk4.cpp

n-body.cpp
```

---

# File Descriptions

## common.hpp

This file contains constants and common includes used throughout the project.

### Libraries

* `<cmath>` – Mathematical functions
* `<iostream>` – Console input/output
* `<fstream>` – Reading and writing files
* `<string>` – String handling
* `<vector>` – Dynamic arrays

### Unit System

The simulation uses astronomical units rather than SI units.

| Quantity | Unit                   |
| -------- | ---------------------- |
| Distance | Astronomical Unit (AU) |
| Mass     | Solar Mass             |
| Time     | Year                   |

Using this system,

[
G = 4\pi^2
]

which simplifies many orbital calculations.

### Constants

* `PI`
* `G`
* `AU`
* `SOLAR_MASS`
* `SOLAR_YEAR`

---

## particle.hpp

This file defines the fundamental data structures used throughout the simulation.

### Vector3D

Represents a three-dimensional vector.

Supported operations include:

* Addition
* Subtraction
* Scalar multiplication
* Scalar division
* Dot product
* Cross product
* Vector magnitude
* Vector normalization

The class is used to represent

* Position
* Velocity
* Acceleration

throughout the simulation.

---

### Particle

Represents a single gravitating body.

Each particle stores:

* Name (optional identifier)
* Mass
* Position
* Velocity
* Acceleration

The acceleration vector is updated every timestep by the selected force calculation algorithm.

---

## force_calculator.hpp

Defines the abstract interface for all force calculation algorithms.

Every force calculator must implement:

* `computeAccelerations()`
* `name()`

The base class also stores the gravitational softening parameter.

This allows different algorithms to be used interchangeably through runtime polymorphism.

Current implementation:

* Naive pairwise algorithm

Future implementations:

* Barnes-Hut Tree
* Fast Multipole Method (optional)

---

## force_calculator_naive

Implements the direct pairwise gravitational interaction.

Algorithm:

1. Reset every particle's acceleration.
2. Iterate over every unique particle pair.
3. Compute the softened gravitational interaction.
4. Update the accelerations of both particles simultaneously using Newton's Third Law.

Time complexity:

* **O(N²)**

Memory complexity:

* **O(1)** (excluding particle storage)

---

## integrator.hpp

Defines the abstract interface for all numerical integration methods.

Every integrator implements:

* `step()`
* `name()`

Each integrator receives

* the particle system,
* a force calculator,
* and the timestep,

allowing any force calculation algorithm to be combined with any integration scheme.

---

## Euler Integrator

Implements the explicit Euler method.

Characteristics:

* First-order accurate
* Simple implementation
* Poor long-term energy conservation
* Mainly included for comparison purposes

---

## Euler-Cromer Integrator

Updates velocity before position.

Compared to Euler:

* Better stability
* Improved energy behaviour
* Still first-order accurate

---

## Leapfrog Integrator

Implements the Kick-Drift-Kick (KDK) Leapfrog scheme.

Characteristics:

* Second-order accurate
* Symplectic
* Excellent long-term energy conservation
* Well suited for gravitational simulations

---

## Velocity Verlet Integrator

Implements the Velocity Verlet method.

Characteristics:

* Second-order accurate
* Symplectic
* Conserves energy well over long simulations
* Widely used in molecular dynamics and astrophysical simulations

---

## RK4 Integrator

Implements the classical fourth-order Runge-Kutta method.

Characteristics:

* Fourth-order accurate
* Requires four force evaluations per timestep
* Very accurate for short integrations
* Not symplectic, therefore exhibits long-term energy drift

---

## n-body.cpp

The main driver of the simulation.

Responsibilities include:

* Creating the particle system
* Initializing particle masses, positions and velocities
* Selecting the force calculation algorithm
* Selecting the numerical integrator
* Executing the simulation loop
* Writing simulation results to output files

---

# Object-Oriented Design

The project is designed around two independent abstractions.

```
ForceCalculator
│
├── NaiveForceCalculator
└── BarnesHutForceCalculator (remaining)

TimeIntegrator
│
├── EulerIntegrator
├── EulerCromerIntegrator
├── LeapfrogIntegrator
├── VelocityVerletIntegrator
└── RK4Integrator
```

This design allows force calculation algorithms and numerical integration methods to be combined independently without modifying the remainder of the simulation.

---

