#pragma once

#include "force_calculator.hpp"

class TimeIntegrator {
    public:
        virtual ~TimeIntegrator() = default;

        virtual void step(
            std::vector<Particle>& particles,
            ForceCalculator& calculator,
            double dt
        ) = 0;

        virtual std::string name() const = 0;
};