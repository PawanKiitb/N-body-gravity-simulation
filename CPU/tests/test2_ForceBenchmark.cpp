#include "../include/tests_common.hpp"

#include <random>
#include <algorithm>

int main()
{
    std::mt19937 rng(42);
    for(int times = 0; times < 5; times++)
    {
        std::cout
        << "\n=================================================================\n";
        std::cout
        << "Run : "
        << times + 1
        << "\n";
        std::cout
        << "===================================================================\n";

    std::uniform_real_distribution<double> position(-100.0, 100.0);
    std::uniform_real_distribution<double> mass(0.5, 2.0);

    std::vector<int> particleCounts =
    {
        100,
        500,
        1000,
        5000,
        10000,
    };
    
    for(int N : particleCounts)
    {
        std::cout
        << "\n=========================================\n";
        std::cout
        << "Number of Particles : "
        << N
        << "\n";
        std::cout
        << "=========================================\n";

        std::vector<Particle> particles;
        particles.reserve(N);

        for(int i = 0; i < N; i++)
        {
            particles.emplace_back(
                mass(rng),
                Vector3D(
                    position(rng),
                    position(rng),
                    position(rng)
                ),
                Vector3D(0.0,0.0,0.0),
                "Particle " + std::to_string(i)
            );
        }        

        std::vector<Particle> naiveParticles = particles;
        std::vector<Particle> bhParticles = particles;
        std::vector<Particle> bhOptimizedParticles = particles;

        ForceCalculatorNaive naive;
        ForceCalculatorBarnesHut bh(0.0,0.5);
        ForceCalculatorBarnesHutOptimized bh2(0.0,0.5, 8);

        for(int i=0;i<5;i++) naive.computeAccelerations(naiveParticles);

        for(int i=0;i<5;i++) bh.computeAccelerations(bhParticles);

        for(int i=0;i<5;i++) bh2.computeAccelerations(bhOptimizedParticles);

        naive.getTimeStatistics().reset();
        bh.getTimeStatistics().reset();
        bh2.getTimeStatistics().reset();

        // Compute accelerations
        for(int i = 0; i<100; i++) naive.computeAccelerations(naiveParticles);
        for(int i = 0; i<100; i++) bh.computeAccelerations(bhParticles);
        for(int i = 0; i<100; i++) bh2.computeAccelerations(bhOptimizedParticles);

        // Compute errors
        double averageAbsoluteError = 0.0;
        double maximumAbsoluteError = 0.0;

        double averageRelativeError = 0.0;
        double maximumRelativeError = 0.0;

        for(size_t i = 0; i < naiveParticles.size(); i++)
        {
            Vector3D diff =
                bhParticles[i].acceleration -
                naiveParticles[i].acceleration;

            double absError = diff.norm();

            averageAbsoluteError += absError;
            maximumAbsoluteError =
                std::max(maximumAbsoluteError, absError);

            double reference =
                naiveParticles[i].acceleration.norm();

            if(reference > 1e-12)
            {
                double relError = absError / reference;

                averageRelativeError += relError;
                maximumRelativeError =
                    std::max(maximumRelativeError, relError);
            }
        }

        averageAbsoluteError /= N;
        averageRelativeError /= N;

        double averageAbsoluteErrorOptimized = 0.0;
        double maximumAbsoluteErrorOptimized = 0.0;

        double averageRelativeErrorOptimized = 0.0;
        double maximumRelativeErrorOptimized = 0.0;

        for(size_t i = 0; i < naiveParticles.size(); i++)
        {
            Vector3D diff =
                bhOptimizedParticles[i].acceleration -
                naiveParticles[i].acceleration;

            double absError = diff.norm();

            averageAbsoluteErrorOptimized += absError;
            maximumAbsoluteErrorOptimized =
                std::max(maximumAbsoluteErrorOptimized, absError);

            double reference =
                naiveParticles[i].acceleration.norm();

            if(reference > 1e-12)
            {
                double relError = absError / reference;

                averageRelativeErrorOptimized += relError;
                maximumRelativeErrorOptimized =
                    std::max(maximumRelativeErrorOptimized, relError);
            }
        }

        averageAbsoluteErrorOptimized /= N;
        averageRelativeErrorOptimized /= N;

        // Timing statistics
        std::cout << "\nNaive\n";
        std::cout
        << "Force Time : "
        << naive.getTimeStatistics().total_time
        << " ms, "
        << "Average    : "
        << naive.getTimeStatistics().average_time()
        << " ms\n";

        std::cout << "\nBarnes-Hut\n";

        std::cout
        << "Total Time : "
        << bh.getTimeStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh.getTimeStatistics().average_time()
        << " ms\n";

        std::cout
        << "Tree Build : "
        << bh.getBuildingTreeStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh.getBuildingTreeStatistics().average_time()
        << " ms\n";

        std::cout
        << "Traversal  : "
        << bh.getTraversalStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh.getTraversalStatistics().average_time()
        << " ms\n";

        std::cout
        << "Destroy    : "
        << bh.getDestroyingTreeStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh.getDestroyingTreeStatistics().average_time()
        << " ms\n";

        std::cout
        << "\nNaive Speedup : "
        << naive.getTimeStatistics().total_time /
           bh.getTimeStatistics().total_time
        << "x\n";

        std::cout
        << "\nAverage Absolute Error : "
        << averageAbsoluteError
        << '\n';

        std::cout
        << "Maximum Absolute Error : "
        << maximumAbsoluteError
        << '\n';

        std::cout
        << "Average Relative Error : "
        << averageRelativeError
        << '\n';

        std::cout
        << "Maximum Relative Error : "
        << maximumRelativeError
        << '\n';

        std::cout << "\nBarnes-Hut Optimized\n";

        std::cout
        << "Total Time : "
        << bh2.getTimeStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh2.getTimeStatistics().average_time()
        << " ms\n";

        std::cout
        << "Tree Build : "
        << bh2.getBuildingTreeStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh2.getBuildingTreeStatistics().average_time()
        << " ms\n";

        std::cout
        << "Traversal  : "
        << bh2.getTraversalStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh2.getTraversalStatistics().average_time()
        << " ms\n";

        std::cout
        << "Destroy    : "
        << bh2.getDestroyingTreeStatistics().total_time
        << " ms, "
        << "Average    : "
        << bh2.getDestroyingTreeStatistics().average_time()
        << " ms\n";

        std::cout
        << "\nNaive Speedup : "
        << naive.getTimeStatistics().total_time /
           bh2.getTimeStatistics().total_time
        << "x\n";

        std::cout
        << "\nBarnes-Hut Speedup : "
        << bh.getTimeStatistics().total_time /
           bh2.getTimeStatistics().total_time
        << "x\n";

        std::cout
        << "\nAverage Absolute Error : "
        << averageAbsoluteErrorOptimized
        << '\n';

        std::cout
        << "Maximum Absolute Error : "
        << maximumAbsoluteErrorOptimized
        << '\n';

        std::cout
        << "Average Relative Error : "
        << averageRelativeErrorOptimized
        << '\n';

        std::cout
        << "Maximum Relative Error : "
        << maximumRelativeErrorOptimized
        << '\n';
    }

    std::cout << "===================================================================\n";
    std::cout << "Run " << times + 1 << " completed.\n";
    std::cout << "===================================================================\n";
}

    return 0;
}