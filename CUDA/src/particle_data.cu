#include "../include/particle_data.hpp"
#include <cstring>
#include <iostream>

// Allocate device memory for particle data and related arrays
void allocateParticles(ParticleArrays &particles, size_t num_particles)
{
    particles.num_particles = num_particles;
    size_t bytes = num_particles * sizeof(double);
    size_t bytes_u32 = num_particles * sizeof(uint32_t);
    size_t bytes_int = num_particles * sizeof(int);

    // Positions and velocities
    cudaMalloc(&particles.pos_x, bytes);
    cudaMalloc(&particles.pos_y, bytes);
    cudaMalloc(&particles.pos_z, bytes);
    cudaMalloc(&particles.vel_x, bytes);
    cudaMalloc(&particles.vel_y, bytes);
    cudaMalloc(&particles.vel_z, bytes);
    cudaMalloc(&particles.acc_x, bytes);
    cudaMalloc(&particles.acc_y, bytes);
    cudaMalloc(&particles.acc_z, bytes);
    cudaMalloc(&particles.old_acc_x, bytes);
    cudaMalloc(&particles.old_acc_y, bytes);
    cudaMalloc(&particles.old_acc_z, bytes);
    cudaMalloc(&particles.mass, bytes);

    // Morton code and indices (for Barnes-Hut)
    cudaMalloc(&particles.morton_codes, bytes_u32);
    cudaMalloc(&particles.morton_codes_sorted, bytes_u32);
    cudaMalloc(&particles.particle_indices, bytes_int);
    cudaMalloc(&particles.particle_indices_sorted, bytes_int);

    // Sorted arrays for positions, mass, and accelerations
    cudaMalloc(&particles.pos_x_sorted, bytes);
    cudaMalloc(&particles.pos_y_sorted, bytes);
    cudaMalloc(&particles.pos_z_sorted, bytes);
    cudaMalloc(&particles.mass_sorted, bytes);
    cudaMalloc(&particles.acc_x_sorted, bytes);
    cudaMalloc(&particles.acc_y_sorted, bytes);
    cudaMalloc(&particles.acc_z_sorted, bytes);
}

// Free device memory
void freeParticles(ParticleArrays &particles)
{
    // Free all allocated arrays
    cudaFree(particles.pos_x);
    cudaFree(particles.pos_y);
    cudaFree(particles.pos_z);

    cudaFree(particles.vel_x);
    cudaFree(particles.vel_y);
    cudaFree(particles.vel_z);

    cudaFree(particles.acc_x);
    cudaFree(particles.acc_y);
    cudaFree(particles.acc_z);

    cudaFree(particles.old_acc_x);
    cudaFree(particles.old_acc_y);
    cudaFree(particles.old_acc_z);

    cudaFree(particles.mass);

    cudaFree(particles.morton_codes);
    cudaFree(particles.morton_codes_sorted);
    cudaFree(particles.particle_indices);
    cudaFree(particles.particle_indices_sorted);

    cudaFree(particles.pos_x_sorted);
    cudaFree(particles.pos_y_sorted);
    cudaFree(particles.pos_z_sorted);
    cudaFree(particles.mass_sorted);
    cudaFree(particles.acc_x_sorted);
    cudaFree(particles.acc_y_sorted);
    cudaFree(particles.acc_z_sorted);

    // Reset particle count and pointers
    particles.num_particles = 0;
    particles.pos_x = particles.pos_y = particles.pos_z = nullptr;
    particles.vel_x = particles.vel_y = particles.vel_z = nullptr;
    particles.acc_x = particles.acc_y = particles.acc_z = nullptr;
    particles.old_acc_x = particles.old_acc_y = particles.old_acc_z = nullptr;
    particles.mass = nullptr;
    particles.morton_codes = particles.morton_codes_sorted = nullptr;
    particles.particle_indices = particles.particle_indices_sorted = nullptr;
    particles.pos_x_sorted = particles.pos_y_sorted = particles.pos_z_sorted = nullptr;
    particles.mass_sorted = nullptr;
    particles.acc_x_sorted = particles.acc_y_sorted = particles.acc_z_sorted = nullptr;
}

// Copy data from host arrays to device arrays
void copyParticlesToDevice(
    ParticleArrays &particles,
    const double *pos_x, const double *pos_y, const double *pos_z,
    const double *vel_x, const double *vel_y, const double *vel_z,
    const double *mass)
{
    size_t bytes = particles.num_particles * sizeof(double);
    // Copy positions and velocities
    cudaMemcpy(particles.pos_x, pos_x, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.pos_y, pos_y, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.pos_z, pos_z, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.vel_x, vel_x, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.vel_y, vel_y, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.vel_z, vel_z, bytes, cudaMemcpyHostToDevice);
    // Copy masses
    cudaMemcpy(particles.mass, mass, bytes, cudaMemcpyHostToDevice);

    // Initialize accelerations and old accelerations to zero on device
    cudaMemset(particles.acc_x, 0, bytes);
    cudaMemset(particles.acc_y, 0, bytes);
    cudaMemset(particles.acc_z, 0, bytes);
    cudaMemset(particles.old_acc_x, 0, bytes);
    cudaMemset(particles.old_acc_y, 0, bytes);
    cudaMemset(particles.old_acc_z, 0, bytes);
}

// Copy results from device back to host arrays
void copyParticlesToHost(
    const ParticleArrays &particles,
    double *pos_x, double *pos_y, double *pos_z,
    double *vel_x, double *vel_y, double *vel_z,
    double *acc_x, double *acc_y, double *acc_z,
    double *mass)
{
    size_t bytes = particles.num_particles * sizeof(double);
    // Copy positions and velocities back
    cudaMemcpy(pos_x, particles.pos_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(pos_y, particles.pos_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(pos_z, particles.pos_z, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_x, particles.vel_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_y, particles.vel_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_z, particles.vel_z, bytes, cudaMemcpyDeviceToHost);
    // Copy accelerations
    cudaMemcpy(acc_x, particles.acc_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(acc_y, particles.acc_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(acc_z, particles.acc_z, bytes, cudaMemcpyDeviceToHost);
    // Copy masses back
    cudaMemcpy(mass, particles.mass, bytes, cudaMemcpyDeviceToHost);
}
