#pragma once

#include "integrator.hpp"

class EulerIntegrator : public TimeIntegrator {
    protected:
        std::string integrator_name = "Euler Integrator";
    public: 
        void step(
            std::vector<Particle>& particles,
            ForceCalculator& calculator,
            double dt
        ) override;

        std::string name() const override;
};