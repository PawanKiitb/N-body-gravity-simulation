#pragma once

#include "../include/integrator_velocity_verlet.hpp"

void VelocityVerletIntegrator::step(std::vector<Particle>& particles, ForceCalculator& calculator, double dt) {
    //  update the position based on current velocity and acceleration
    calculator.computeAccelerations(particles);
    for (Particle& p : particles) {
        p.position += p.velocity * dt + 0.5 * p.acceleration * dt * dt;
    }

    // update the vlocity based on the new accelerations
    // store the old accelerations
    old_accelerations.reserve(particles.size());
    for (const Particle& p : particles) {
        old_accelerations.push_back(p.acceleration);
    }
    // compute the new accelerations
    calculator.computeAccelerations(particles);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].velocity += (old_accelerations[i] + particles[i].acceleration) * (0.5 * dt);
    }
}

std::string VelocityVerletIntegrator::name() const {
    return "Velocity-Verlet";
}