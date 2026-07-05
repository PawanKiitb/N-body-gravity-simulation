#pragma once

#include "../include/common.hpp"
#include "../include/particle_data.hpp"
#include "../include/force_calculator_naive.hpp"
#include "../include/integrator_euler.hpp"
#include "../include/integrator_euler_cromer.hpp"
#include "../include/integrator_leapfrog.hpp"
#include "../include/integrator_velocity_verlet.hpp"
#include "../include/integrator_rk4.hpp"

#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>

inline void copyStateToHost(
    const ParticleArrays& particles,
    std::vector<double>& pos_x,
    std::vector<double>& pos_y,
    std::vector<double>& pos_z,
    std::vector<double>& vel_x,
    std::vector<double>& vel_y,
    std::vector<double>& vel_z,
    std::vector<double>& acc_x,
    std::vector<double>& acc_y,
    std::vector<double>& acc_z,
    std::vector<double>& mass
)
{
    size_t n = particles.num_particles;

    pos_x.resize(n);
    pos_y.resize(n);
    pos_z.resize(n);

    vel_x.resize(n);
    vel_y.resize(n);
    vel_z.resize(n);

    acc_x.resize(n);
    acc_y.resize(n);
    acc_z.resize(n);

    mass.resize(n);

    copyParticlesToHost(
        particles,
        pos_x.data(),
        pos_y.data(),
        pos_z.data(),
        vel_x.data(),
        vel_y.data(),
        vel_z.data(),
        acc_x.data(),
        acc_y.data(),
        acc_z.data(),
        mass.data()
    );
}

inline double kineticEnergy(
    const std::vector<double>& vel_x,
    const std::vector<double>& vel_y,
    const std::vector<double>& vel_z,
    const std::vector<double>& mass
)
{
    double ke = 0.0;

    for(size_t i=0;i<mass.size();i++)
    {
        ke += 0.5 * mass[i] *
            (vel_x[i]*vel_x[i]
           + vel_y[i]*vel_y[i]
           + vel_z[i]*vel_z[i]);
    }

    return ke;
}

inline double potentialEnergy(
    const std::vector<double>& pos_x,
    const std::vector<double>& pos_y,
    const std::vector<double>& pos_z,
    const std::vector<double>& mass
)
{
    double pe = 0.0;

    size_t n = mass.size();

    for(size_t i=0;i<n;i++)
    {
        for(size_t j=i+1;j<n;j++)
        {
            double dx = pos_x[j]-pos_x[i];
            double dy = pos_y[j]-pos_y[i];
            double dz = pos_z[j]-pos_z[i];

            double r = std::sqrt(dx*dx+dy*dy+dz*dz);

            pe -= constant::G * mass[i] * mass[j] / r;
        }
    }

    return pe;
}

inline double totalEnergy(
    const std::vector<double>& pos_x,
    const std::vector<double>& pos_y,
    const std::vector<double>& pos_z,
    const std::vector<double>& vel_x,
    const std::vector<double>& vel_y,
    const std::vector<double>& vel_z,
    const std::vector<double>& mass
)
{
    return kineticEnergy(
                vel_x,
                vel_y,
                vel_z,
                mass)
         + potentialEnergy(
                pos_x,
                pos_y,
                pos_z,
                mass);
}

inline void angularMomentum(
    const std::vector<double>& pos_x,
    const std::vector<double>& pos_y,
    const std::vector<double>& pos_z,
    const std::vector<double>& vel_x,
    const std::vector<double>& vel_y,
    const std::vector<double>& vel_z,
    const std::vector<double>& mass,
    double& Lx,
    double& Ly,
    double& Lz
)
{
    Lx = Ly = Lz = 0.0;

    for(size_t i=0;i<mass.size();i++)
    {
        Lx += mass[i] * (pos_y[i]*vel_z[i]-pos_z[i]*vel_y[i]);
        Ly += mass[i] * (pos_z[i]*vel_x[i]-pos_x[i]*vel_z[i]);
        Lz += mass[i] * (pos_x[i]*vel_y[i]-pos_y[i]*vel_x[i]);
    }
}

inline void centerOfMass(
    const std::vector<double>& pos_x,
    const std::vector<double>& pos_y,
    const std::vector<double>& pos_z,
    const std::vector<double>& mass,
    double& cx,
    double& cy,
    double& cz
)
{
    cx = cy = cz = 0.0;

    double totalMass = 0.0;

    for(size_t i=0;i<mass.size();i++)
    {
        cx += mass[i] * pos_x[i];
        cy += mass[i] * pos_y[i];
        cz += mass[i] * pos_z[i];

        totalMass += mass[i];
    }

    cx /= totalMass;
    cy /= totalMass;
    cz /= totalMass;
}

inline double earthSunDistance(
    const std::vector<double>& pos_x,
    const std::vector<double>& pos_y,
    const std::vector<double>& pos_z
)
{
    double dx = pos_x[1]-pos_x[0];
    double dy = pos_y[1]-pos_y[0];
    double dz = pos_z[1]-pos_z[0];

    return std::sqrt(dx*dx+dy*dy+dz*dz);
}