#pragma once

#include "integrator.hpp"

class EulerCromerIntegrator : public TimeIntegrator {
    protected:
        std::string integrator_name = "Euler Integrator";
        cudaEvent_t startEvent;
        cudaEvent_t stopEvent;

        double total_compute_time_ms;
        bool is_first_step = true; // Flag to check if it's the first step
    public: 
        EulerCromerIntegrator() : total_compute_time_ms(0.0) {
            cudaEventCreate(&startEvent);
            cudaEventCreate(&stopEvent);
        }
        ~EulerCromerIntegrator() {
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