#pragma once

#include "integrator.hpp"

class EulerIntegrator : public TimeIntegrator {
    protected:
        std::string integrator_name = "Euler Integrator";
        cudaEvent_t startEvent;
        cudaEvent_t stopEvent;

        double total_compute_time_ms;
    public: 
        EulerIntegrator() : total_compute_time_ms(0.0) {
            cudaEventCreate(&startEvent);
            cudaEventCreate(&stopEvent);
        }
        ~EulerIntegrator() {
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