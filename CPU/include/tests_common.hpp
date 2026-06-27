#pragma once
#include "common.hpp"
#include "particle.hpp"
#include "force_calculator_naive.hpp"
#include "force_calculator_barnes_hut.hpp"
#include "force_calculator_barnes_hut_optimized.hpp"
#include "integrator_euler.hpp"
#include "integrator_velocity_verlet.hpp"
#include "integrator_rk4.hpp"
#include "integrator_euler-cromer.hpp"
#include "integrator_leapfrog.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <memory>

template<typename Integrator, typename Calculator>
std::vector<Particle> runSimulation(
    std::vector<Particle> particles,
    Integrator& integrator,
    Calculator& calculator,
    double dt,
    int steps)
{
    for(int i=0;i<steps;i++)
        integrator.step(particles,calculator,dt);

    return particles;
}

double kineticEnergy(const std::vector<Particle>& particles)
{
    double KE = 0.0;
    for (const Particle& p : particles) {
        KE += 0.5 * p.mass * p.velocity.dot(p.velocity);
    }
    return KE;
}

double potentialEnergy(const std::vector<Particle>& particles) {
    double PE = 0.0;
    size_t n = particles.size();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            Vector3D r = particles[j].position - particles[i].position;
            double distance = r.norm();
            PE -= constant::G * particles[i].mass * particles[j].mass / distance;
        }
    }
    return PE;
}

double totalEnergy(const std::vector<Particle>& particles) {
    return kineticEnergy(particles) + potentialEnergy(particles);
}

Vector3D angularMomentum(const std::vector<Particle>& particles) {
    Vector3D L;
    for (const Particle& p : particles) {
        L += p.position.cross(p.velocity)*p.mass;
    }
    return L;
}

Vector3D centerOfMass(const std::vector<Particle>& particles) {
    Vector3D com;
    double totalMass = 0.0;
    for (const Particle& p : particles) {
        com += p.position * p.mass;
        totalMass += p.mass;
    }
    return com / totalMass;
}

double positionError(const std::vector<Particle>& a,
                     const std::vector<Particle>& b)
{
    double err = 0.0;

    for(size_t i=0;i<a.size();i++)
    {
        err += (a[i].position-b[i].position).norm();
    }

    return err/a.size();
}

double velocityError(const std::vector<Particle>& a,
                     const std::vector<Particle>& b)
{
    double err = 0.0;

    for(size_t i=0;i<a.size();i++)
    {
        err += (a[i].velocity-b[i].velocity).norm();
    }

    return err/a.size();
}