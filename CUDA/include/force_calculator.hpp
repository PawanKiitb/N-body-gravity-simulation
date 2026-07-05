#pragma once

#include <common.hpp>
#include <particle_data.hpp>

/*
Base class for computing acceleration on each particle using CUDA
The base class will have these methods:
- Constructors: Default and one with softening parameter
- computeAccelerations: Return the acceleration on each particle(with Softening)
- name: Return the name of the method used for computing acceleration
*/

class ForceCalculator {
    protected:
        double softening_squared;
    public:
        // constructor
        ForceCalculator() : softening_squared(0.0) {}
        ForceCalculator(double softening) : softening_squared(softening*softening) {}

        // destructor
        virtual ~ForceCalculator() = default;

        virtual void computeAccelerations(
            ParticleArrays& particles
        ) = 0;

        virtual std::string name() const = 0;
};