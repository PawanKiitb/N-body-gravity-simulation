#include "../include/tests_common.hpp"

#include <random>
#include <algorithm>
#include <iomanip>

struct Result
{
    int leafCapacity;

    double averageTime = 0.0;

    double averageRelativeError = 0.0;
    double maximumRelativeError = 0.0;
};

int main()
{
    constexpr int NUM_EXPERIMENTS = 5;
    constexpr int NUM_ITERATIONS  = 100;

    std::vector<int> particleCounts =
    {
        100,
        500,
        1000,
        5000,
        10000
    };

    std::vector<int> leafCapacities =
    {
        1,
        2,
        4,
        8,
        16,
        32,
        64
    };

    for(int N : particleCounts)
    {
        std::cout
        << "\n=============================================================\n";
        std::cout
        << "Particle Count : "
        << N
        << "\n";
        std::cout
        << "=============================================================\n";

        //----------------------------------------------------------
        // One accumulated result for each leaf capacity
        //----------------------------------------------------------

        std::vector<Result> results;

        for(int leaf : leafCapacities)
        {
            Result r;
            r.leafCapacity = leaf;
            results.push_back(r);
        }

        //----------------------------------------------------------
        // Repeat experiment several times
        //----------------------------------------------------------

        for(int experiment = 0;
            experiment < NUM_EXPERIMENTS;
            experiment++)
        {
            std::mt19937 rng(42 + experiment);

            std::uniform_real_distribution<double> position(-100.0,100.0);
            std::uniform_real_distribution<double> mass(0.5,2.0);

            //------------------------------------------------------
            // Generate particles
            //------------------------------------------------------

            std::vector<Particle> particles;
            particles.reserve(N);

            for(int i=0;i<N;i++)
            {
                particles.emplace_back(
                    mass(rng),
                    Vector3D(
                        position(rng),
                        position(rng),
                        position(rng)
                    ),
                    Vector3D(0,0,0),
                    "Particle"
                );
            }

            //------------------------------------------------------
            // Naive reference
            //------------------------------------------------------

            std::vector<Particle> reference = particles;

            ForceCalculatorNaive naive;

            naive.computeAccelerations(reference);

            //------------------------------------------------------
            // Test each leaf capacity
            //------------------------------------------------------

            for(size_t k=0;k<leafCapacities.size();k++)
            {
                int leaf = leafCapacities[k];

                std::vector<Particle> current = particles;

                ForceCalculatorBarnesHutOptimized bh(
                    0.0,
                    0.5,
                    leaf
                );

                //---------------- Warmup ----------------

                for(int i=0;i<5;i++)
                    bh.computeAccelerations(current);

                bh.getTimeStatistics().reset();

                //---------------- Benchmark ------------

                for(int i=0;i<NUM_ITERATIONS;i++)
                    bh.computeAccelerations(current);

                //---------------- Error ----------------

                double avgRel = 0.0;
                double maxRel = 0.0;

                for(int i=0;i<N;i++)
                {
                    double ref =
                        reference[i].acceleration.norm();

                    if(ref < 1e-12)
                        continue;

                    double rel =
                        (current[i].acceleration -
                         reference[i].acceleration).norm()
                        / ref;

                    avgRel += rel;

                    maxRel =
                        std::max(maxRel, rel);
                }

                avgRel /= N;

                results[k].averageTime +=
                    bh.getTimeStatistics().average_time();

                results[k].averageRelativeError +=
                    avgRel;

                results[k].maximumRelativeError +=
                    maxRel;
            }
        }

        //----------------------------------------------------------
        // Average over experiments
        //----------------------------------------------------------

        for(Result& r : results)
        {
            r.averageTime /= NUM_EXPERIMENTS;
            r.averageRelativeError /= NUM_EXPERIMENTS;
            r.maximumRelativeError /= NUM_EXPERIMENTS;
        }

        //----------------------------------------------------------
        // Sort by runtime
        //----------------------------------------------------------

        auto byTime = results;

        std::sort(
            byTime.begin(),
            byTime.end(),
            [](const Result& a,
               const Result& b)
            {
                return a.averageTime <
                       b.averageTime;
            });

        std::cout
        << "\n================ Fastest =================\n\n";

        std::cout
        << std::left
        << std::setw(8)  << "Leaf"
        << std::setw(15) << "Time(ms)"
        << std::setw(18) << "Avg Rel Error"
        << std::setw(18) << "Max Rel Error"
        << "\n";

        for(const Result& r : byTime)
        {
            std::cout
            << std::setw(8)
            << r.leafCapacity

            << std::setw(15)
            << r.averageTime

            << std::setw(18)
            << r.averageRelativeError

            << std::setw(18)
            << r.maximumRelativeError

            << "\n";
        }

        //----------------------------------------------------------
        // Sort by accuracy
        //----------------------------------------------------------

        auto byAccuracy = results;

        std::sort(
            byAccuracy.begin(),
            byAccuracy.end(),
            [](const Result& a,
               const Result& b)
            {
                return a.averageRelativeError <
                       b.averageRelativeError;
            });

        std::cout
        << "\n============== Most Accurate ==============\n\n";

        std::cout
        << std::left
        << std::setw(8)  << "Leaf"
        << std::setw(18) << "Avg Rel Error"
        << std::setw(18) << "Max Rel Error"
        << std::setw(15) << "Time(ms)"
        << "\n";

        for(const Result& r : byAccuracy)
        {
            std::cout
            << std::setw(8)
            << r.leafCapacity

            << std::setw(18)
            << r.averageRelativeError

            << std::setw(18)
            << r.maximumRelativeError

            << std::setw(15)
            << r.averageTime

            << "\n";
        }
    }

    return 0;
}