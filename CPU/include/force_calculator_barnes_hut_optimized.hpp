#pragma once

#include "force_calculator.hpp"

class ForceCalculatorBarnesHutOptimized : public ForceCalculator {
    protected:
        double theta_sqr;
        int rootIndex;
        std::vector<Node> nodePool;
        size_t nextFreeNode = 0;
        std::vector<Particle>* particles_ = nullptr;
        int leaf_capacity;
    public:
        TimingStatistics total = {
            .name = "Barnes-Hut Force Calculation",
            .total_time = 0.0,
            .calls = 0
        };
        TimingStatistics buildingTree = {
            .name = "Barnes-Hut Tree Building",
            .total_time = 0.0,
            .calls = 0
        };
        TimingStatistics traversal = {
            .name = "Barnes-Hut Tree Traversal",
            .total_time = 0.0,
            .calls = 0
        };
        TimingStatistics destroyingTree = {
            .name = "Barnes-Hut Tree Destruction",
            .total_time = 0.0,
            .calls = 0
        };

        double epsilon = 1e-9;
        // constructor
        ForceCalculatorBarnesHutOptimized() : ForceCalculator(), theta_sqr(0.5*0.5), rootIndex(-1) {
            nodePool.resize(constant::MAX_NODES);
        }
        ForceCalculatorBarnesHutOptimized(double softening, double theta, int leaf_capacity) : ForceCalculator(softening), theta_sqr(theta*theta), rootIndex(-1), leaf_capacity(leaf_capacity) {
            nodePool.resize(constant::MAX_NODES);
        }

        // destructor
        ~ForceCalculatorBarnesHutOptimized() = default;

        // main method to compute accelerations
        void computeAccelerations(
            std::vector<Particle>& particles
        ) override;

        // tree methods
        void buildTree(std::vector<Particle>& particles);
        void insertParticle(int nodeIndex, int particleIndex);

        // tree helper methods
        int allocateNode();
        bool subdivideNode(int nodeIndex);
        void destroyTree();

        // Force calculation methods
        void computeAccelerationForParticle(int nodeIndex, int particleIndex);

        // other methods
        std::string name() const override;
        TimingStatistics& getTimeStatistics() const override;
        TimingStatistics& getBuildingTreeStatistics() const {
            return const_cast<TimingStatistics&>(buildingTree);
        }
        TimingStatistics& getTraversalStatistics() const {
            return const_cast<TimingStatistics&>(traversal);
        }
        TimingStatistics& getDestroyingTreeStatistics() const {
            return const_cast<TimingStatistics&>(destroyingTree);
        }
};