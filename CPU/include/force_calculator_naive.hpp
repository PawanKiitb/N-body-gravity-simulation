#pragma once

#include "force_calculator.hpp"

class ForceCalculatorNaive : public ForceCalculator {
public:
    TimingStatistics naive_timing = {
        .name = "Naive Force Calculation",
        .total_time = 0.0,
        .calls = 0
    };
    using ForceCalculator::ForceCalculator;

    ~ForceCalculatorNaive() override = default;

    void computeAccelerations(
        std::vector<Particle>& particles
    ) override;

    std::string name() const override;

    TimingStatistics& getTimeStatistics() const override;
};