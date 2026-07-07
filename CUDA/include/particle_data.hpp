#pragma once

#include <common.hpp>
#include <cstdint>

struct ParticleArrays
{
    int num_particles;

    double *pos_x;
    double *pos_y;
    double *pos_z;

    double *vel_x;
    double *vel_y;
    double *vel_z;

    double *acc_x;
    double *acc_y;
    double *acc_z;

    double *old_acc_x;
    double *old_acc_y;
    double *old_acc_z;

    double *mass;

    uint32_t *morton_codes;
    uint32_t *morton_codes_sorted;

    int *particle_indices;
    int *particle_indices_sorted;
};

void allocateParticles(
    ParticleArrays &particles,
    size_t num_particles);

void freeParticles(
    ParticleArrays &particles);

void copyParticlesToDevice(
    ParticleArrays &particles,
    const double *pos_x,
    const double *pos_y,
    const double *pos_z,
    const double *vel_x,
    const double *vel_y,
    const double *vel_z,
    const double *mass);

void copyParticlesToHost(
    const ParticleArrays &particles,
    double *pos_x,
    double *pos_y,
    double *pos_z,
    double *vel_x,
    double *vel_y,
    double *vel_z,
    double *acc_x,
    double *acc_y,
    double *acc_z,
    double *mass);