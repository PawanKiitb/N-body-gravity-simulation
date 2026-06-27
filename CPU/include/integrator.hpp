#pragma once

#include "force_calculator.hpp"

class TimeIntegrator {
    protected:
        std::string integrator_name;
    public:
        TimingStatistics timing_statistics = {
            .name = integrator_name,
            .total_time = 0.0,
            .calls = 0
        };
        virtual ~TimeIntegrator() = default;

        virtual void step(
            std::vector<Particle>& particles,
            ForceCalculator& calculator,
            double dt
        ) = 0;

        virtual std::string name() const = 0;
        virtual const TimingStatistics& getTimeStatistics() const {
            return timing_statistics;
        }
};