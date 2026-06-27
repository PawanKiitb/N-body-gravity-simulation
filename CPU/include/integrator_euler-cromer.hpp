#pragma once

#include "integrator.hpp"

class EulerCromerIntegrator : public TimeIntegrator {
    protected:
        std::string integrator_name = "Euler-Cromer Integrator";
    public:
        void step(
            std::vector<Particle>& particles,
            ForceCalculator& calculator,
            double dt
        ) override;
        std::string name() const override;
};