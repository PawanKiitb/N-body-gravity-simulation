#include "../include/tests_common.hpp"

int main() {
    // create a csv file to write the output, which will store x, y and z
    std::ofstream output_file("orbit.csv");
    output_file << "earth_x,earth_y,earth_z,sun_x,sun_y,sun_z,kinetic_energy,potential_energy,total_energy,angular_momentum_x,angular_momentum_y,angular_momentum_z,center_of_mass_x,center_of_mass_y,center_of_mass_z\n";

    // Let's try to simulate Earth and Sun system
    Particle sun(1, Vector3D(0, 0, 0), Vector3D(0, 0, 0), "Sun");
    double earth_mass = 3.003e-6; // Mass of Earth in solar masses;
    Particle earth(earth_mass, Vector3D(1, 0, 0), Vector3D(0, 2*constant::PI, 0),"Earth");
    std::vector<Particle> particles = {sun, earth};
    ForceCalculatorBarnesHut force_calculator(0.0, 0.3);
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