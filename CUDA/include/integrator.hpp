#pragma once

#include <particle_data.hpp>
#include <force_calculator.hpp>

class Integrator {
    public:
        virtual ~Integrator() = default;
        virtual void step(
            ForceCalculator& calculator,
            ParticleArrays& particles,
            double dt
        ) = 0;

    virtual std::string name() const = 0;
};
