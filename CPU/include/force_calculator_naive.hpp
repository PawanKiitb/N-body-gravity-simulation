#pragma once

#include "force_calculator.hpp"

class NaiveForceCalculator : public ForceCalculator {
public:
    using ForceCalculator::ForceCalculator;

    ~NaiveForceCalculator() override = default;

    void computeAccelerations(
        std::vector<Particle>& particles
    ) override;

    std::string name() const override;
};