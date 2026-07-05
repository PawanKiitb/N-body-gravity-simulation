#include "../include/particle_data.hpp"

#include <cstring>
#include <iostream>

void allocateParticles(ParticleArrays& particles, size_t num_particles)
{
    particles.num_particles = num_particles;

    cudaMalloc(&particles.pos_x, num_particles * sizeof(double));
    cudaMalloc(&particles.pos_y, num_particles * sizeof(double));
    cudaMalloc(&particles.pos_z, num_particles * sizeof(double));

    cudaMalloc(&particles.vel_x, num_particles * sizeof(double));
    cudaMalloc(&particles.vel_y, num_particles * sizeof(double));
    cudaMalloc(&particles.vel_z, num_particles * sizeof(double));

    cudaMalloc(&particles.acc_x, num_particles * sizeof(double));
    cudaMalloc(&particles.acc_y, num_particles * sizeof(double));
    cudaMalloc(&particles.acc_z, num_particles * sizeof(double));

    cudaMalloc(&particles.old_acc_x, num_particles * sizeof(double));
    cudaMalloc(&particles.old_acc_y, num_particles * sizeof(double));
    cudaMalloc(&particles.old_acc_z, num_particles * sizeof(double));

    cudaMalloc(&particles.mass, num_particles * sizeof(double));
}

void freeParticles(ParticleArrays& particles)
{
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

    particles.num_particles = 0;

    particles.pos_x = nullptr;
    particles.pos_y = nullptr;
    particles.pos_z = nullptr;

    particles.vel_x = nullptr;
    particles.vel_y = nullptr;
    particles.vel_z = nullptr;

    particles.acc_x = nullptr;
    particles.acc_y = nullptr;
    particles.acc_z = nullptr;

    particles.old_acc_x = nullptr;
    particles.old_acc_y = nullptr;
    particles.old_acc_z = nullptr;

    particles.mass = nullptr;
}

void copyParticlesToDevice(
    ParticleArrays& particles,
    const double* pos_x,
    const double* pos_y,
    const double* pos_z,
    const double* vel_x,
    const double* vel_y,
    const double* vel_z,
    const double* mass
)
{
    size_t bytes = particles.num_particles * sizeof(double);

    cudaMemcpy(particles.pos_x, pos_x, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.pos_y, pos_y, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.pos_z, pos_z, bytes, cudaMemcpyHostToDevice);

    cudaMemcpy(particles.vel_x, vel_x, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.vel_y, vel_y, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(particles.vel_z, vel_z, bytes, cudaMemcpyHostToDevice);

    cudaMemcpy(particles.mass, mass, bytes, cudaMemcpyHostToDevice);

    cudaMemset(particles.acc_x, 0, bytes);
    cudaMemset(particles.acc_y, 0, bytes);
    cudaMemset(particles.acc_z, 0, bytes);

    cudaMemset(particles.old_acc_x, 0, bytes);
    cudaMemset(particles.old_acc_y, 0, bytes);
    cudaMemset(particles.old_acc_z, 0, bytes);
}

void copyParticlesToHost(
    const ParticleArrays& particles,
    double* pos_x,
    double* pos_y,
    double* pos_z,
    double* vel_x,
    double* vel_y,
    double* vel_z,
    double* acc_x,
    double* acc_y,
    double* acc_z,
    double* mass
)
{
    size_t bytes = particles.num_particles * sizeof(double);

    cudaMemcpy(pos_x, particles.pos_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(pos_y, particles.pos_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(pos_z, particles.pos_z, bytes, cudaMemcpyDeviceToHost);

    cudaMemcpy(vel_x, particles.vel_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_y, particles.vel_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(vel_z, particles.vel_z, bytes, cudaMemcpyDeviceToHost);

    cudaMemcpy(acc_x, particles.acc_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(acc_y, particles.acc_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(acc_z, particles.acc_z, bytes, cudaMemcpyDeviceToHost);

    cudaMemcpy(mass, particles.mass, bytes, cudaMemcpyDeviceToHost);
}