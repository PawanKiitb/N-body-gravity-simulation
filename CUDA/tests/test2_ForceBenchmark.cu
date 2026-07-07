#include "../src/test_common.cpp"

#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

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
        10000
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

            ParticleArrays naiveParticles;
            ParticleArrays bhParticles;

            allocateParticles(naiveParticles, N);
            allocateParticles(bhParticles, N);

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
                h_mass[i] = massDist(rng);
            }

            copyParticlesToDevice(
                naiveParticles,
                h_pos_x.data(),
                h_pos_y.data(),
                h_pos_z.data(),
                h_vel_x.data(),
                h_vel_y.data(),
                h_vel_z.data(),
                h_mass.data()
            );

            copyParticlesToDevice(
                bhParticles,
                h_pos_x.data(),
                h_pos_y.data(),
                h_pos_z.data(),
                h_vel_x.data(),
                h_vel_y.data(),
                h_vel_z.data(),
                h_mass.data()
            );

            ForceCalculatorNaive naive(0.0);
            BarnesHutForceCalculator bh(0.0, 0.5);

            // Warm-up
            for(int i = 0; i < 5; i++)
            {
                naive.computeAccelerations(naiveParticles);
                bh.computeAccelerations(bhParticles);
            }

            naive.resetTimeStatistics();
            bh.resetTimeStatistics();

            for(int i = 0; i < 100; i++)
                naive.computeAccelerations(naiveParticles);

            for(int i = 0; i < 100; i++)
                bh.computeAccelerations(bhParticles);

            std::vector<double> naive_pos_x;
            std::vector<double> naive_pos_y;
            std::vector<double> naive_pos_z;

            std::vector<double> naive_vel_x;
            std::vector<double> naive_vel_y;
            std::vector<double> naive_vel_z;

            std::vector<double> naive_acc_x;
            std::vector<double> naive_acc_y;
            std::vector<double> naive_acc_z;

            std::vector<double> naive_mass;

            std::vector<double> bh_pos_x;
            std::vector<double> bh_pos_y;
            std::vector<double> bh_pos_z;

            std::vector<double> bh_vel_x;
            std::vector<double> bh_vel_y;
            std::vector<double> bh_vel_z;

            std::vector<double> bh_acc_x;
            std::vector<double> bh_acc_y;
            std::vector<double> bh_acc_z;

            std::vector<double> bh_mass;

            copyStateToHost(
                naiveParticles,
                naive_pos_x,
                naive_pos_y,
                naive_pos_z,
                naive_vel_x,
                naive_vel_y,
                naive_vel_z,
                naive_acc_x,
                naive_acc_y,
                naive_acc_z,
                naive_mass
            );

            copyStateToHost(
                bhParticles,
                bh_pos_x,
                bh_pos_y,
                bh_pos_z,
                bh_vel_x,
                bh_vel_y,
                bh_vel_z,
                bh_acc_x,
                bh_acc_y,
                bh_acc_z,
                bh_mass
            );

            double avgAbs = 0.0;
            double maxAbs = 0.0;

            double avgRel = 0.0;
            double maxRel = 0.0;

            for(int i = 0; i < N; i++)
            {
                double dx = bh_acc_x[i] - naive_acc_x[i];
                double dy = bh_acc_y[i] - naive_acc_y[i];
                double dz = bh_acc_z[i] - naive_acc_z[i];

                double absErr = std::sqrt(dx*dx + dy*dy + dz*dz);

                avgAbs += absErr;
                maxAbs = std::max(maxAbs, absErr);

                double ref =
                    std::sqrt(
                        naive_acc_x[i]*naive_acc_x[i] +
                        naive_acc_y[i]*naive_acc_y[i] +
                        naive_acc_z[i]*naive_acc_z[i]);

                if(ref > 1e-12)
                {
                    double rel = absErr / ref;

                    avgRel += rel;
                    maxRel = std::max(maxRel, rel);
                }
            }

            avgAbs /= N;
            avgRel /= N;

            std::cout << "\nNaive\n";

            std::cout
                << "Total Time : "
                << naive.getTotalComputeTime()
                << " ms\n";

            std::cout << "\nBarnes-Hut\n";

            std::cout
                << "Total Time : "
                << bh.getTotalComputeTime()
                << " ms\n";

            std::cout
                << "Tree Build : "
                << bh.getBuildTreeTime()
                << " ms\n";

            std::cout
                << "Traversal  : "
                << bh.getTraversalTime()
                << " ms\n";

            std::cout
                << "Destroy    : "
                << bh.getDestroyTreeTime()
                << " ms\n";

            std::cout
                << "\nSpeedup : "
                << naive.getTotalComputeTime() /
                   bh.getTotalComputeTime()
                << "x\n";

            std::cout
                << "\nAverage Absolute Error : "
                << avgAbs << '\n';

            std::cout
                << "Maximum Absolute Error : "
                << maxAbs << '\n';

            std::cout
                << "Average Relative Error : "
                << avgRel << '\n';

            std::cout
                << "Maximum Relative Error : "
                << maxRel << '\n';

            freeParticles(naiveParticles);
            freeParticles(bhParticles);
        }
    }

    return 0;
}