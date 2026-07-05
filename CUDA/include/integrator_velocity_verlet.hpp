#pragma once

#include "integrator.hpp"

class VelocityVerletIntegrator : public TimeIntegrator {
    protected:
        cudaEvent_t startEvent, stopEvent;
        double total_compute_time_ms;
        bool is_first_step = true; 
    public: 
        VelocityVerletIntegrator() : total_compute_time_ms(0.0) {
            cudaEventCreate(&startEvent);
            cudaEventCreate(&stopEvent);
        }
        ~VelocityVerletIntegrator() {
            cudaEventDestroy(startEvent);
            cudaEventDestroy(stopEvent);
        }
        void step(
            ForceCalculator& calculator,
            ParticleArrays& particles,
            double dt
        ) override;
        std::string name() const override;
};