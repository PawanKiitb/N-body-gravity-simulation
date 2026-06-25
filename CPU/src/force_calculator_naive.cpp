#pragma once

#include "../include/force_calculator_naive.hpp"

void NaiveForceCalculator::computeAccelerations(std::vector<Particle>& particles) {
    for (Particle& p : particles) {
        p.acceleration = {0.0, 0.0, 0.0};
    }
    size_t n = particles.size();
    for(size_t i = 0; i<n; i++) {
        for(size_t j = i+1; j<n; j++) {
            Vector3D r_ij = particles[j].position - particles[i].position;
            double distance_squared = r_ij.dot(r_ij) + softening_squared;
            double distance = std::sqrt(distance_squared);
            double distance_cubed = distance * distance * distance;
            particles[i].acceleration = particles[i].acceleration + r_ij * (constant::G * particles[j].mass / distance_cubed);
            particles[j].acceleration = particles[j].acceleration - r_ij * (constant::G * particles[i].mass / distance_cubed);
        }
    }
}

std::string NaiveForceCalculator::name() const {
    return "Naive";
}