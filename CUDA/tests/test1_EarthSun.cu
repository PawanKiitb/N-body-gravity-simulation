#include "../src/test_common.cpp"

int main()
{
    constexpr int NUM_PARTICLES = 2;

    ParticleArrays particles;
    allocateParticles(particles, NUM_PARTICLES);

    //-----------------------------
    // Host-side initialization
    //-----------------------------

    double pos_x[NUM_PARTICLES] = {0.0, 1.0};
    double pos_y[NUM_PARTICLES] = {0.0, 0.0};
    double pos_z[NUM_PARTICLES] = {0.0, 0.0};

    double vel_x[NUM_PARTICLES] = {0.0, 0.0};
    double vel_y[NUM_PARTICLES] = {0.0, 2.0 * constant::PI};
    double vel_z[NUM_PARTICLES] = {0.0, 0.0};

    double mass[NUM_PARTICLES] = {
        1.0,
        3.003e-6
    };

    copyParticlesToDevice(
        particles,
        pos_x,
        pos_y,
        pos_z,
        vel_x,
        vel_y,
        vel_z,
        mass
    );

    //--------------------------------
    // Choose calculator/integrator
    //--------------------------------

    ForceCalculatorNaive calculator(0.0);

    VelocityVerletIntegrator integrator;
    // EulerIntegrator integrator;
    // EulerCromerIntegrator integrator;
    // LeapfrogIntegrator integrator;
    // RK4Integrator integrator;

    //--------------------------------
    // CSV
    //--------------------------------

    std::ofstream file("orbit.csv");

    file
    << "earth_x,"
    << "earth_y,"
    << "earth_z,"
    << "sun_x,"
    << "sun_y,"
    << "sun_z,"
    << "earth_distance,"
    << "kinetic_energy,"
    << "potential_energy,"
    << "total_energy,"
    << "angular_momentum_x,"
    << "angular_momentum_y,"
    << "angular_momentum_z,"
    << "center_of_mass_x,"
    << "center_of_mass_y,"
    << "center_of_mass_z\n";

    //--------------------------------
    // Temporary host arrays
    //--------------------------------

    std::vector<double> h_pos_x;
    std::vector<double> h_pos_y;
    std::vector<double> h_pos_z;

    std::vector<double> h_vel_x;
    std::vector<double> h_vel_y;
    std::vector<double> h_vel_z;

    std::vector<double> h_acc_x;
    std::vector<double> h_acc_y;
    std::vector<double> h_acc_z;

    std::vector<double> h_mass;

    //--------------------------------
    // Simulation
    //--------------------------------

    double dt = 1.0 / 365.0;
    int steps = static_cast<int>(1.0 / dt);

    for(int step = 0; step <= steps; step++)
    {
        copyStateToHost(
            particles,
            h_pos_x,
            h_pos_y,
            h_pos_z,
            h_vel_x,
            h_vel_y,
            h_vel_z,
            h_acc_x,
            h_acc_y,
            h_acc_z,
            h_mass
        );

        double ke = kineticEnergy(
            h_vel_x,
            h_vel_y,
            h_vel_z,
            h_mass
        );

        double pe = potentialEnergy(
            h_pos_x,
            h_pos_y,
            h_pos_z,
            h_mass
        );

        double total = ke + pe;

        double Lx, Ly, Lz;

        angularMomentum(
            h_pos_x,
            h_pos_y,
            h_pos_z,
            h_vel_x,
            h_vel_y,
            h_vel_z,
            h_mass,
            Lx,
            Ly,
            Lz
        );

        double com_x, com_y, com_z;

        centerOfMass(
            h_pos_x,
            h_pos_y,
            h_pos_z,
            h_mass,
            com_x,
            com_y,
            com_z
        );

        file
            << std::fixed
            << std::setprecision(15)

            << h_pos_x[1] << ","
            << h_pos_y[1] << ","
            << h_pos_z[1] << ","

            << h_pos_x[0] << ","
            << h_pos_y[0] << ","
            << h_pos_z[0] << ","

            << earthSunDistance(
                h_pos_x,
                h_pos_y,
                h_pos_z
            ) << ","

            << ke << ","
            << pe << ","
            << total << ","

            << Lx << ","
            << Ly << ","
            << Lz << ","

            << com_x << ","
            << com_y << ","
            << com_z << "\n";

        integrator.step(
            calculator,
            particles,
            dt
        );
    }

    file.close();

    freeParticles(particles);

    std::cout << "Simulation complete.\n";
    std::cout << "Output written to orbit.csv\n";

    return 0;
}