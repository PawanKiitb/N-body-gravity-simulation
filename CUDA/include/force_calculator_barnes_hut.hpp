#pragma once

#include "force_calculator.hpp"

// Axis-aligned bounding box of the system
struct Bounds {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double size;
};

// Structure holding the Barnes–Hut tree node data (all on device)
struct TreeArrays {
    double *center_x, *center_y, *center_z;
    double *size, *size_sqr;
    double *mass;
    double *com_x, *com_y, *com_z;
    int    *first_child;    // index of first child node (-1 if leaf)
    int    *parent;         // index of parent node
    int    *particle_count; // number of particles in this node (leaf nodes)
    int    *particle_start; // start index in the sorted-particle array
    int    *node_depth;     // depth of this node in the tree
};

struct BarnesHutForceCalculator : public ForceCalculator {
    double softening_squared;
    double theta_sqr;
    int leaf_capacity;

    // Timing events
    cudaEvent_t startEvent, stopEvent;
    double total_compute_time_ms;
    double build_tree_time_ms;
    double destroy_tree_time_ms;
    double traverse_tree_time_ms;

    // Bounding box (on device)
    Bounds* d_bounding_cube;
    void* bounding_box_temp_storage = nullptr;
    size_t bounding_box_temp_storage_bytes = 0;

    // Radix sort temporary storage for Morton codes
    void* radix_sort_temp_storage = nullptr;
    size_t radix_sort_temp_storage_bytes = 0;

    // Tree arrays and management
    TreeArrays tree;
    int* nextFreeNode;  // pointer to next free node index on device

public:
    BarnesHutForceCalculator(double softening, double theta, int leaf_capacity = constant::LEAF_CAPACITY)
        : softening_squared(softening * softening), theta_sqr(theta*theta),
          leaf_capacity(leaf_capacity), total_compute_time_ms(0.0),
          build_tree_time_ms(0.0), destroy_tree_time_ms(0.0), traverse_tree_time_ms(0.0)
    {
        cudaEventCreate(&startEvent);
        cudaEventCreate(&stopEvent);

        // Allocate memory for bounding box (device)
        cudaMalloc((void**)&d_bounding_cube, sizeof(Bounds));

        // Allocate memory for tree arrays (device)
        cudaMalloc((void**)&tree.center_x, sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.center_y, sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.center_z, sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.size,     sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.size_sqr, sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.mass,     sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.com_x,    sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.com_y,    sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.com_z,    sizeof(double) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.first_child,    sizeof(int) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.parent,         sizeof(int) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.particle_count, sizeof(int) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.particle_start, sizeof(int) * constant::MAX_NODES);
        cudaMalloc((void**)&tree.node_depth,     sizeof(int) * constant::MAX_NODES);

        // Allocate memory for next-free-node counter
        cudaMalloc((void**)&nextFreeNode, sizeof(int));
    }

    ~BarnesHutForceCalculator() {
        // Destroy CUDA events
        cudaEventDestroy(startEvent);
        cudaEventDestroy(stopEvent);
        // Free all CUDA memory allocated in constructor
        cudaFree(d_bounding_cube);
        cudaFree(tree.center_x);     cudaFree(tree.center_y);     cudaFree(tree.center_z);
        cudaFree(tree.size);         cudaFree(tree.size_sqr);
        cudaFree(tree.mass);
        cudaFree(tree.com_x);        cudaFree(tree.com_y);        cudaFree(tree.com_z);
        cudaFree(tree.first_child);
        cudaFree(tree.parent);
        cudaFree(tree.particle_count); cudaFree(tree.particle_start);
        cudaFree(tree.node_depth);
        cudaFree(nextFreeNode);
        // Free temporary storages if allocated
        if (bounding_box_temp_storage) cudaFree(bounding_box_temp_storage);
        if (radix_sort_temp_storage)  cudaFree(radix_sort_temp_storage);
    }

    // Compute accelerations for all particles (implements ForceCalculator)
    void computeAccelerations(ParticleArrays& particles) override;

    // Tree-building helpers
    void buildTree(ParticleArrays& particles);
    void computeBoundingBox(ParticleArrays& particles);
    void radixSort(ParticleArrays& particles);
    void buildOctree(const uint32_t* sorted_morton_codes, int num_particles, TreeArrays& tree, int* d_next_free_node);
    void computeMassCOM(ParticleArrays& particles);
    void destroyTree();

    std::string name() const override { return "Barnes-Hut"; }
    double getTotalComputeTime() const { return total_compute_time_ms; }
    double getBuildTreeTime() const { return build_tree_time_ms; }
    double getDestroyTreeTime() const { return destroy_tree_time_ms; }
    double getTraversalTime() const { return traverse_tree_time_ms; }
    void resetTotalComputeTime() { total_compute_time_ms = 0.0; build_tree_time_ms = 0.0; destroy_tree_time_ms = 0.0; traverse_tree_time_ms = 0.0; }
};
