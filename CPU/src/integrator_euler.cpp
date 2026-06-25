#include "../include/integrator_euler.hpp"

void EulerIntegrator::step(std::vector<Particle>& particles, ForceCalculator& calculator, double dt) {
    // compute the accelerations using the provided force calculator
    calculator.computeAccelerations(particles);
    
    // update the position
    for (Particle& p : particles) {
        p.position = p.position + p.velocity * dt;
    }
    // update the velocity
    for (Particle& p : particles) {
        p.velocity = p.velocity + p.acceleration * dt;
    }
}

std::string EulerIntegrator::name() const {
    return "Euler";
}