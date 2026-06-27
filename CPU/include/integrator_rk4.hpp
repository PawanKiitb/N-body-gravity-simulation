#pragma once

#include "integrator.hpp"

class RK4Integrator : public TimeIntegrator {
    private:
        std::vector<Particle> K2State;
        std::vector<Particle> K3State;
        std::vector<Particle> K4State;
    protected:
        std::string integrator_name = "RK4 Integrator";
    public: 
        void step(
            std::vector<Particle>& particles,
            ForceCalculator& calculator,
            double dt
        ) override;
        std::string name() const override;
};