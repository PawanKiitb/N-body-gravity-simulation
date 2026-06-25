#pragma once

#include "../include/integrator_rk4.hpp"

void RK4Integrator::step(std::vector<Particle>& particles, ForceCalculator& calculator, double dt) {
    // Copy the particles state
    size_t n = particles.size();
    K2State.resize(n);
    K3State.resize(n);
    K4State.resize(n);
    K2State = particles;
    K3State = particles;
    K4State = particles;

    // Step 1: Compute K1
    calculator.computeAccelerations(particles);

    // step 2: Compute K2
    for (size_t i = 0; i < n; i++) {
        K2State[i].position = particles[i].position + particles[i].velocity * (dt / 2.0);
        K2State[i].velocity = particles[i].velocity + particles[i].acceleration * (dt / 2.0);
    }
    calculator.computeAccelerations(K2State);

    // step 3: Compute K3
    for (size_t i = 0; i < n; i++) {
        K3State[i].position = particles[i].position + K2State[i].velocity * (dt / 2.0);
        K3State[i].velocity = particles[i].velocity + K2State[i].acceleration * (dt / 2.0);
    }
    calculator.computeAccelerations(K3State);

    // step 4: Compute K4
    for (size_t i = 0; i < n; i++) {
        K4State[i].position = particles[i].position + K3State[i].velocity * dt;
        K4State[i].velocity = particles[i].velocity + K3State[i].acceleration * dt;
    }
    calculator.computeAccelerations(K4State);

    // Update the particles state
    for (size_t i = 0; i < n; i++) {
        particles[i].position += (particles[i].velocity + 2.0 * K2State[i].velocity + 2.0 * K3State[i].velocity + K4State[i].velocity) * (dt / 6.0);
        particles[i].velocity += (particles[i].acceleration + 2.0 * K2State[i].acceleration + 2.0 * K3State[i].acceleration + K4State[i].acceleration) * (dt / 6.0);
    }
}

std::string RK4Integrator::name() const {
    return "RK4";
}