#pragma once

#include "integrator.hpp"

class LeapfrogIntegrator : public TimeIntegrator {
    protected:
        std::string integrator_name = "Leapfrog Integrator";
    public: 
    void step(
        std::vector<Particle>& particles,
        ForceCalculator& calculator,
        double dt
    ) override;

    std::string name() const override;
};