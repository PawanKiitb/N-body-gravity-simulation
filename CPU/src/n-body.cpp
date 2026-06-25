#include "../include/common.hpp"
#include "../include/particle.hpp"
#include "../include/force_calculator_naive.hpp"
#include "../include/integrator_euler.hpp"
#include "../include/integrator_velocity_verlet.hpp"
#include "../include/integrator_rk4.hpp"
#include "../include/integrator_euler-cromer.hpp"
#include "../include/integrator_leapfrog.hpp"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>

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

int main() {
    // create a csv file to write the output, which will store x, y and z
    std::ofstream output_file("orbit.csv");
    output_file << "earth_x,earth_y,earth_z,sun_x,sun_y,sun_z,kinetic_energy,potential_energy,total_energy,angular_momentum_x,angular_momentum_y,angular_momentum_z,center_of_mass_x,center_of_mass_y,center_of_mass_z\n";

    // Let's try to simulate Earth and Sun system
    Particle sun(1, Vector3D(0, 0, 0), Vector3D(0, 0, 0), "Sun");
    double earth_mass = 3.003e-6; // Mass of Earth in solar masses;
    Particle earth(earth_mass, Vector3D(1, 0, 0), Vector3D(0, 2*constant::PI, 0),"Earth");
    std::vector<Particle> particles = {sun, earth};
    NaiveForceCalculator force_calculator(0.0);
    VelocityVerletIntegrator integrator;
    
    Vector3D L = angularMomentum(particles);
    Vector3D com = centerOfMass(particles);

    output_file
    << particles[1].position.x << ','
    << particles[1].position.y << ','
    << particles[1].position.z << ','
    << particles[0].position.x << ','
    << particles[0].position.y << ','
    << particles[0].position.z << ','
    << kineticEnergy(particles) << ','
    << potentialEnergy(particles) << ','
    << totalEnergy(particles) << ','
    << L.x << ','
    << L.y << ','
    << L.z << ','
    << com.x << ','
    << com.y << ','
    << com.z << '\n';
    
    double dt = 1.0/365.0;
    int steps = static_cast<int>(1.0 / dt);
    for (int step = 0; step <= steps; step++) {
        integrator.step(particles, force_calculator, dt);
        // write to csv file
        L = angularMomentum(particles);
        com = centerOfMass(particles);
        output_file << std::fixed << std::setprecision(15) 
                    << particles[1].position.x << ","
                    << particles[1].position.y << ","
                    << particles[1].position.z << ","
                    << particles[0].position.x << ","
                    << particles[0].position.y << ","
                    << particles[0].position.z << ","
                    << kineticEnergy(particles) << ","
                    << potentialEnergy(particles) << ","
                    << totalEnergy(particles) << ","
                    << L.x << ","
                    << L.y << ","
                    << L.z << ","
                    << com.x << ","
                    << com.y << ","
                    << com.z << "\n";
        
    }
    particles[0].info();
    std::cout << "\n";
    particles[1].info();
    output_file.close();
    return 0;
}