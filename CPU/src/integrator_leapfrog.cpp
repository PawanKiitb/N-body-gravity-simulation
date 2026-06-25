#pragma once

#include "../include/integrator_leapfrog.hpp"

void LeapfrogIntegrator::step(std::vector<Particle>& particles, ForceCalculator& calculator, double dt) {
    calculator.computeAccelerations(particles);
    double half_dt = dt / 2.0;
    // update half step velocity
    for (Particle& p : particles) {
        p.velocity += p.acceleration * half_dt;
    }

    // update position based on these half step velocities
    for (Particle& p : particles) {
        p.position += p.velocity * dt;
    }

    // update velocity based on new accelerations
    calculator.computeAccelerations(particles);
    for (Particle& p : particles) {
        p.velocity += p.acceleration * half_dt;
    }
}

std::string LeapfrogIntegrator::name() const {
    return "Leapfrog";
}