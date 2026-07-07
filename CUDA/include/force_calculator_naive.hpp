#pragma once

#include "force_calculator.hpp"

class ForceCalculatorNaive : public ForceCalculator {
    protected:
        cudaEvent_t startEvent;
        cudaEvent_t stopEvent;

        double total_compute_time_ms;
    public:
        using ForceCalculator::ForceCalculator;
        ForceCalculatorNaive(double softening) : ForceCalculator(softening), total_compute_time_ms(0.0) {
            cudaEventCreate(&startEvent);
            cudaEventCreate(&stopEvent);
        }
        ~ForceCalculatorNaive() override = default;

        void computeAccelerations(
            ParticleArrays& particles
        ) override;

        std::string name() const override;
        // funtion to return the total compute time in milliseconds
        double getTotalComputeTime() const {
            return total_compute_time_ms;
        }
};