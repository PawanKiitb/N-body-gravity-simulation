#include "../include/tests_common.hpp"

int main() {
    double dt = 1.0/365.0;
    int steps = 365;

    // Reference
    ForceCalculatorNaive naiveReference;
    VelocityVerletIntegrator verletReference;

    std::vector<Particle> particles = {
        Particle(1.0, Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 0.0, 0.0), "Sun"),
        Particle(3.003e-6, Vector3D(1.0, 0.0, 0.0), Vector3D(0.0, 2*M_PI, 0.0), "Earth")
    };

    auto reference = 
    runSimulation(
            particles,
            verletReference,
            naiveReference,
            dt,
            steps);
    
    std::vector<std::unique_ptr<TimeIntegrator>> integrators;

    integrators.emplace_back(std::make_unique<EulerIntegrator>());
    integrators.emplace_back(std::make_unique<EulerCromerIntegrator>());
    integrators.emplace_back(std::make_unique<LeapfrogIntegrator>());
    integrators.emplace_back(std::make_unique<VelocityVerletIntegrator>());
    integrators.emplace_back(std::make_unique<RK4Integrator>());

    for(auto& integrator : integrators)
{
    std::cout
    << "\n=====================================\n";

    std::cout
    << integrator->name()
    << "\n";

    ForceCalculatorNaive naive;

    auto naiveParticles =
        runSimulation(
            particles,
            *integrator,
            naive,
            dt,
            steps);

    std::cout
    << "\nNaive\n";

    const auto& n =
        naive.getTimeStatistics();

    std::cout
    << "Force time : "
    << n.total_time
    << " ms\n";

    std::cout
    << "Average    : "
    << n.average_time()
    << " ms\n";

    std::cout
    << "Position error : "
    << positionError(
        naiveParticles,
        reference)
    << '\n';

    std::cout
    << "Velocity error : "
    << velocityError(
        naiveParticles,
        reference)
    << '\n';

    ForceCalculatorBarnesHut bh(0.0,0.3);

    auto bhParticles =
        runSimulation(
            particles,
            *integrator,
            bh,
            dt,
            steps);

    std::cout
    << "\nBarnes-Hut\n";

    std::cout
    << "Force time : "
    << bh.getTimeStatistics().total_time
    << " ms\n";

    std::cout
    << "Tree build : "
    << bh.getBuildingTreeStatistics().total_time
    << " ms\n";

    std::cout
    << "Traversal  : "
    << bh.getTraversalStatistics().total_time
    << " ms\n";

    std::cout
    << "Tree destroy : "
    << bh.getDestroyingTreeStatistics().total_time
    << " ms\n";

    std::cout
    << "Position error : "
    << positionError(
        bhParticles,
        reference)
    << '\n';

    std::cout
    << "Velocity error : "
    << velocityError(
        bhParticles,
        reference)
    << '\n';

    std::cout
    << "=====================================\n";
}


}