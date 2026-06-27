#include "../include/force_calculator_naive.hpp"

void ForceCalculatorNaive::computeAccelerations(std::vector<Particle>& particles) {
    for (Particle& p : particles) {
        p.acceleration = {0.0, 0.0, 0.0};
    }
    size_t n = particles.size();
    auto start = Clock::now();
    // for(size_t i = 0; i<n; i++) {
    //     for(size_t j = i+1; j<n; j++) {
    //         Vector3D r_ij = particles[j].position - particles[i].position;
    //         double distance_squared = r_ij.dot(r_ij) + softening_squared;
    //         double inv_distance = 1/std::sqrt(distance_squared);
    //         double inv_distance_cubed = inv_distance * inv_distance * inv_distance;
    //         particles[i].acceleration = particles[i].acceleration + r_ij * (constant::G * particles[j].mass * inv_distance_cubed);
    //         particles[j].acceleration = particles[j].acceleration - r_ij * (constant::G * particles[i].mass * inv_distance_cubed);
    //     }
    // }

    #pragma omp parallel for schedule(static)
    for(size_t i = 0; i<n; i++) {
        for(size_t j = 0; j<n; j++) {
            if(i==j) continue;
            Vector3D r_ij = particles[j].position - particles[i].position;
            double distance_squared = r_ij.dot(r_ij) + softening_squared;
            double inv_distance = 1/std::sqrt(distance_squared);
            double inv_distance_cubed = inv_distance * inv_distance * inv_distance;
            particles[i].acceleration += r_ij * (constant::G * particles[j].mass * inv_distance_cubed);
        }
    }
    auto end = Clock::now();
    double elapsed_time = std::chrono::duration<double, std::milli>(end - start).count();
    naive_timing.calls++;
    naive_timing.total_time += elapsed_time;
}

std::string ForceCalculatorNaive::name() const {
    return "Naive";
}

TimingStatistics& ForceCalculatorNaive::getTimeStatistics() const {
    return const_cast<TimingStatistics&>(naive_timing);
}