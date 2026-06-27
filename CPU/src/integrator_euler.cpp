#include "../include/integrator_euler.hpp"

void EulerIntegrator::step(std::vector<Particle>& particles, ForceCalculator& calculator, double dt) {
    auto start = Clock::now();
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
    auto end = Clock::now();
    double elapsed_time = std::chrono::duration<double, std::milli>(end - start).count();
    timing_statistics.calls++;
    timing_statistics.total_time += elapsed_time;
}

std::string EulerIntegrator::name() const {
    return "Euler";
}