#pragma once

#include "force_calculator.hpp"

class ForceCalculatorBarnesHut : public ForceCalculator {
    protected:
        double theta_sqr;
        Node_* root;
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
        ForceCalculatorBarnesHut() : ForceCalculator(), theta_sqr(0.5*0.5), root(nullptr) {}
        ForceCalculatorBarnesHut(double softening, double theta) : ForceCalculator(softening), theta_sqr(theta*theta), root(nullptr) {}

        // destructor
        ~ForceCalculatorBarnesHut() = default;
        void computeAccelerations(
            std::vector<Particle>& particles
        ) override;
        // tree methods
        void buildTree(std::vector<Particle>& particles);
        void insertParticle(Node_* node, Particle* particle);
        int getChildIndex(Node_* node, const Vector3D& particle_position);
        bool subdivideNode(Node_* node);
        void destroyTree(Node_* node);
        // Force calculation methods
        void computeAccelerationForParticle(Node_* node, Particle* particle);
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