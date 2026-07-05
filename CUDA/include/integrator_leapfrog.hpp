#pragma once

#include "integrator.hpp"

class LeapfrogIntegrator : public TimeIntegrator {
    protected:
        cudaEvent_t startEvent, stopEvent;
        double total_compute_time_ms;
        bool is_first_step = true; 
    public: 
    LeapfrogIntegrator() : total_compute_time_ms(0.0) {
        cudaEventCreate(&startEvent);
        cudaEventCreate(&stopEvent);
    }
    ~LeapfrogIntegrator() {
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