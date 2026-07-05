#pragma once

#inlcude <common.hpp>
#include <particle_data.hpp>
#include <force_calculator.hpp>

class TimeIntegrator {
    public:
        virtual ~TimeIntegrator() = default;
        virtual void step(
            ForceCalculator& calculator,
            ParticleArrays& particles,
            double dt
        ) = 0;

    virtual std::string name() const = 0;
};
