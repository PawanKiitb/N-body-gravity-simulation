#include "../include/integrator_euler-cromer.hpp"

void EulerCromerIntegrator::step(std::vector<Particle>& particles, ForceCalculator& calculator, double dt) {
    // compute the accelerations using the provided force calculator
    calculator.computeAccelerations(particles);

    // update the velocity
    for (Particle& p : particles) {
        p.velocity = p.velocity + p.acceleration * dt;
    }
    // update the position
    for (Particle& p : particles) {
        p.position = p.position + p.velocity * dt;
    }
}

std::string EulerCromerIntegrator::name() const {
    return "Euler-Cromer";
}
