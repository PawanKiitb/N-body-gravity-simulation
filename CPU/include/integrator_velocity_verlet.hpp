#pragma once

#include "integrator.hpp"

class VelocityVerletIntegrator : public TimeIntegrator {
    private:
        std::vector<Vector3D> old_accelerations;
    protected:
        std::string integrator_name = "Velocity Verlet Integrator";
    public: 
        void step(
            std::vector<Particle>& particles,
            ForceCalculator& calculator,
            double dt
        ) override;

        std::string name() const override;
};