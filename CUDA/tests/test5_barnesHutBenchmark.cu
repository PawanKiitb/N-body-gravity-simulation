#include "../src/test_common.cpp"

#include <random>
#include <vector>
#include <iostream>

int main()
{
    std::mt19937 rng(42);

    std::uniform_real_distribution<double> position(-100.0, 100.0);
    std::uniform_real_distribution<double> massDist(0.5, 2.0);

    std::vector<int> particleCounts =
    {
        100,
        500,
        1000,
        5000,
        10000,
        50000,
        100000
    };

    for(int run = 0; run < 2; run++)
    {
        std::cout
            << "\n=====================================================\n";
        std::cout
            << "Run " << run + 1 << '\n';
        std::cout
            << "=====================================================\n";

        for(int N : particleCounts)
        {
            std::cout
                << "\n=========================================\n";
            std::cout
                << "Particles : " << N << '\n';
            std::cout
                << "=========================================\n";

            ParticleArrays particles;
            allocateParticles(particles, N);

            std::vector<double> h_pos_x(N);
            std::vector<double> h_pos_y(N);
            std::vector<double> h_pos_z(N);

            std::vector<double> h_vel_x(N, 0.0);
            std::vector<double> h_vel_y(N, 0.0);
            std::vector<double> h_vel_z(N, 0.0);

            std::vector<double> h_mass(N);

            for(int i = 0; i < N; i++)
            {
                h_pos_x[i] = position(rng);
                h_pos_y[i] = position(rng);
                h_pos_z[i] = position(rng);
                h_mass[i]  = massDist(rng);
            }

            copyParticlesToDevice(
                particles,
                h_pos_x.data(),
                h_pos_y.data(),
                h_pos_z.data(),
                h_vel_x.data(),
                h_vel_y.data(),
                h_vel_z.data(),
                h_mass.data()
            );

            BarnesHutForceCalculator bh(0.0, 0.5);

            // Warm-up
            for(int i = 0; i < 5; i++)
                bh.computeAccelerations(particles);

            bh.resetTimeStatistics();

            // Benchmark
            for(int i = 0; i < 100; i++)
                bh.computeAccelerations(particles);

            std::cout << "\nBarnes-Hut\n";

            std::cout
                << "Total Time : "
                << bh.getTotalComputeTime()
                << " ms, Average : "
                << bh.getTotalComputeTime() / 100.0
                << " ms\n";

            std::cout
                << "Tree Build : "
                << bh.getBuildTreeTime()
                << " ms, Average : "
                << bh.getBuildTreeTime() / 100.0
                << " ms\n";

            std::cout
                << "Traversal  : "
                << bh.getTraversalTime()
                << " ms, Average : "
                << bh.getTraversalTime() / 100.0
                << " ms\n";

            std::cout
                << "Destroy    : "
                << bh.getDestroyTreeTime()
                << " ms, Average : "
                << bh.getDestroyTreeTime() / 100.0
                << " ms\n";

            freeParticles(particles);
        }
    }

    return 0;
}