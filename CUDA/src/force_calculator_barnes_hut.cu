#include "../include/force_calculator_barnes_hut.hpp"
#include <iomanip>
#include <fstream>

#define CUDA_CHECK(call)                                                     \
    do {                                                                     \
        cudaError_t err = call;                                              \
        if (err != cudaSuccess) {                                            \
            fprintf(stderr, "CUDA error: %s (err_num=%d)\\n",                \
                    cudaGetErrorString(err), err);                           \
            exit(err);                                                       \
        }                                                                    \
    } while(0)

__global__
void computeAccelerationKernel(
    ParticleArrays particles,
    TreeArrays tree,
    double theta_sqr,
    double softening_squared
)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= particles.num_particles) return;

    float px = particles.pos_x_sorted[tid];
    float py = particles.pos_y_sorted[tid];
    float pz = particles.pos_z_sorted[tid];
    float ax = 0.0f, ay = 0.0f, az = 0.0f;

    constexpr int STACK_SIZE = 128;
    int stack[STACK_SIZE];
    int top = 0;
    stack[top++] = 0;

    while(top > 0) {
        int node = stack[--top];
        if(tree.mass[node] == 0.0f) continue;

        if(tree.first_child[node] == -1) {
            int first = tree.particle_start[node];
            int count = tree.particle_count[node];
            for(int i=0;i<count;i++) {
                int other = first + i;              // direct sorted-space index, no lookup needed
                if(other == tid) continue;
                float dx = particles.pos_x_sorted[other] - px;
                float dy = particles.pos_y_sorted[other] - py;
                float dz = particles.pos_z_sorted[other] - pz;
                float dist2 = dx*dx + dy*dy + dz*dz + static_cast<float>(softening_squared);
                float invDist = rsqrtf(dist2);
                float invDist3 = invDist * invDist * invDist;
                float scale = static_cast<float>(constant::G) * particles.mass_sorted[other] * invDist3;
                ax += dx*scale; ay += dy*scale; az += dz*scale;
            }
            continue;
        }

        int containingChild = 0;
        if(px >= tree.center_x[node]) containingChild |= 4;
        if(py >= tree.center_y[node]) containingChild |= 2;
        if(pz >= tree.center_z[node]) containingChild |= 1;

        int firstChild = tree.first_child[node];

        #pragma unroll
        for(int i=0;i<8;i++) {
            int child = firstChild + i;
            if(tree.mass[child] == 0.0f) continue;

            float dx = tree.com_x[child] - px;
            float dy = tree.com_y[child] - py;
            float dz = tree.com_z[child] - pz;
            float dist2 = dx*dx + dy*dy + dz*dz;

            bool mustDescend = (i == containingChild) ||
                                (tree.size_sqr[child] >= static_cast<float>(theta_sqr) * dist2);

            if (mustDescend && top < STACK_SIZE) {
                stack[top++] = child;
            } else {
                float invDist = rsqrtf(dist2 + static_cast<float>(softening_squared));
                float invDist3 = invDist * invDist * invDist;
                float scale = static_cast<float>(constant::G) * tree.mass[child] * invDist3;
                ax += dx*scale; ay += dy*scale; az += dz*scale;
            }
        }
    }
    particles.acc_x_sorted[tid] = ax;
    particles.acc_y_sorted[tid] = ay;
    particles.acc_z_sorted[tid] = az;
}

__global__
void indexAccelerationKernel(
    ParticleArrays particles,
    const float* acc_x_sorted, const float* acc_y_sorted, const float* acc_z_sorted
)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= particles.num_particles) return;
    int sorted_idx = particles.particle_indices_sorted[tid];
    particles.acc_x[sorted_idx] = static_cast<double>(acc_x_sorted[tid]);
    particles.acc_y[sorted_idx] = static_cast<double>(acc_y_sorted[tid]);
    particles.acc_z[sorted_idx] = static_cast<double>(acc_z_sorted[tid]);
}

// Compute axis-aligned bounding box (min/max) of positions on device
void BarnesHutForceCalculator::computeBoundingBox(ParticleArrays& particles) {
    int n = particles.num_particles;
    if (bounding_box_temp_storage == nullptr) {
        bounding_box_temp_storage_bytes = 0;

        cub::DeviceReduce::Min(
            nullptr,
            bounding_box_temp_storage_bytes,
            particles.pos_x,
            &d_bounding_cube->xmin,
            n
        );

        cudaMalloc(&bounding_box_temp_storage, bounding_box_temp_storage_bytes);
    }

    cub::DeviceReduce::Min(
        bounding_box_temp_storage,
        bounding_box_temp_storage_bytes,
        particles.pos_x,
        &d_bounding_cube->xmin,
        n
    );

    cub::DeviceReduce::Max(
        bounding_box_temp_storage,
        bounding_box_temp_storage_bytes,
        particles.pos_x,
        &d_bounding_cube->xmax,
        n
    );

    // ymin / ymax
    cub::DeviceReduce::Min(
        bounding_box_temp_storage,
        bounding_box_temp_storage_bytes,
        particles.pos_y,
        &d_bounding_cube->ymin,
        n
    );

    cub::DeviceReduce::Max(
        bounding_box_temp_storage,
        bounding_box_temp_storage_bytes,
        particles.pos_y,
        &d_bounding_cube->ymax,
        n
    );

    // zmin / zmax
    cub::DeviceReduce::Min(
        bounding_box_temp_storage,
        bounding_box_temp_storage_bytes,
        particles.pos_z,
        &d_bounding_cube->zmin,
        n
    );

    cub::DeviceReduce::Max(
        bounding_box_temp_storage,
        bounding_box_temp_storage_bytes,
        particles.pos_z,
        &d_bounding_cube->zmax,
        n
    );
}


__global__
void initializeRootKernel(
    Bounds* d_bounding_cube,
    TreeArrays tree,
    int* nextFreeNode
)
{
    // Single thread initializes root node (index 0)
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        float xmin = d_bounding_cube->xmin;
        float ymin = d_bounding_cube->ymin;
        float zmin = d_bounding_cube->zmin;
        float xmax = d_bounding_cube->xmax;
        float ymax = d_bounding_cube->ymax;
        float zmax = d_bounding_cube->zmax;

        float size = fmax(xmax - xmin, fmax(ymax - ymin, zmax - zmin)) + constant::EPSILON;
        d_bounding_cube->size = size;

        tree.center_x[0] = xmin + 0.5 * size;
        tree.center_y[0] = ymin + 0.5 * size;
        tree.center_z[0] = zmin + 0.5 * size;

        tree.size[0] = size;
        tree.size_sqr[0] = size * size;
        tree.mass[0] = 0.0f;
        tree.com_x[0] = 0.0f;
        tree.com_y[0] = 0.0f;
        tree.com_z[0] = 0.0f;
        tree.first_child[0] = -1;
        tree.particle_count[0] = 0;
        tree.particle_start[0] = -1;
        tree.node_depth[0] = 0;
        tree.parent[0] = -1;

        // Initialize next free node index to 1
        *nextFreeNode = 1;
    }
}


__global__
void computeMortonCodesKernel(
    const double* pos_x, const double* pos_y, const double* pos_z,
    const Bounds* bounding_cube,
    uint32_t* morton_codes,
    int* particle_indices,
    int num_particles
)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_particles) return;

    // Compute cube size with epsilon
    float size = bounding_cube->size;

    // Normalize positions to [0,1]
    float nx = (static_cast<float>(pos_x[tid]) - bounding_cube->xmin) / size;
    float ny = (static_cast<float>(pos_y[tid]) - bounding_cube->ymin) / size;
    float nz = (static_cast<float>(pos_z[tid]) - bounding_cube->zmin) / size;

    // Convert to 10-bit integer coordinates
    uint32_t ix = min(max(static_cast<int>(nx * 1024.0), 0), 1023);
    uint32_t iy = min(max(static_cast<int>(ny * 1024.0), 0), 1023);
    uint32_t iz = min(max(static_cast<int>(nz * 1024.0), 0), 1023);

    // Interleave bits to form Morton code
    uint32_t code = 0;
    for(int i = 0; i < 10; i++) {
        code |= ((ix >> i) & 1) << (3*i + 2);
        code |= ((iy >> i) & 1) << (3*i + 1);
        code |= ((iz >> i) & 1) << (3*i);
    }
    morton_codes[tid] = code;
    particle_indices[tid] = tid;
}


void BarnesHutForceCalculator::radixSort(ParticleArrays& particles) {
    if (radix_sort_temp_storage == nullptr) {
        radix_sort_temp_storage_bytes = 0;
        // Query storage size
        cub::DeviceRadixSort::SortPairs(
            nullptr, radix_sort_temp_storage_bytes,
            particles.morton_codes, particles.morton_codes_sorted,
            particles.particle_indices, particles.particle_indices_sorted,
            particles.num_particles
        );
        cudaMalloc(&radix_sort_temp_storage, radix_sort_temp_storage_bytes);
    }
    // Sort the pairs
    cub::DeviceRadixSort::SortPairs(
        radix_sort_temp_storage, radix_sort_temp_storage_bytes,
        particles.morton_codes, particles.morton_codes_sorted,
        particles.particle_indices, particles.particle_indices_sorted,
        particles.num_particles
    );
}

__global__
void sortParticleDataKernel(
    const int* sorted_indices,
    const double* pos_x, const double* pos_y, const double* pos_z,
    const double* mass,
    float* pos_x_sorted, float* pos_y_sorted, float* pos_z_sorted,
    float* mass_sorted,
    int num_particles
)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_particles) return;
    int sorted_idx = sorted_indices[tid];
    pos_x_sorted[tid] = static_cast<float>(pos_x[sorted_idx]);
    pos_y_sorted[tid] = static_cast<float>(pos_y[sorted_idx]);
    pos_z_sorted[tid] = static_cast<float>(pos_z[sorted_idx]);
    mass_sorted[tid]  = static_cast<float>(mass[sorted_idx]);
}

__global__
void buildOctreeLevelKernel(
    const uint32_t* sorted_codes,
    const int* queueFirst, const int* queueCount, const int* queueNodeIdx,
    int queueSize,
    int depth,
    float* center_x, float* center_y, float* center_z,
    float* size, float* size_sqr,
    int* first_child, int* parent,
    int* particle_start, int* particle_count,
    int* d_next_free_node,
    int* queueFirstNext, int* queueCountNext, int* queueNodeIdxNext,
    int* d_next_level_count,
    int leaf_capacity,
    int* node_depth, float* mass,
    float* com_x, float* com_y, float* com_z
) {
    // Your Morton encoder uses 10 bits per axis => 10 tree levels total.
    // Depth 0 uses the topmost octant bits, depth 9 uses the last octant bits.
    constexpr int MORTON_LEVELS = 10;

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= queueSize) return;

    int first    = queueFirst[tid];
    int count    = queueCount[tid];
    int parentId = queueNodeIdx[tid];
    int end      = first + count;

    // If this node should stay a leaf, do not split it.
    if (count <= leaf_capacity || depth >= MORTON_LEVELS - 1) {
        first_child[parentId] = -1;
        particle_start[parentId] = first;
        particle_count[parentId] = count;
        return;
    }

    // Determine which particles fall into each child octant at this depth.
    int childStart[8];
    int childCount[8];

    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        childStart[i] = -1;
        childCount[i] = 0;
    }

    // Extract the 3 bits corresponding to this tree level.
    // With your Morton encoding, depth 0 uses the highest 3 bits,
    // depth 9 uses the lowest 3 bits.
    int mortonShift = 3 * (MORTON_LEVELS - 1 - depth);

    for (int i = first; i < end; ++i) {
        int octant = (sorted_codes[i] >> mortonShift) & 7;
        if (childStart[octant] == -1) {
            childStart[octant] = i;
        }
        childCount[octant]++;
    }

    // Allocate exactly 8 children for this parent.
    int childBase = atomicAdd(d_next_free_node, 8);

    first_child[parentId] = childBase;
    particle_start[parentId] = -1;
    particle_count[parentId] = 0;

    float parentSize = size[parentId];
    float childSize  = parentSize * 0.5;
    float offset     = parentSize * 0.25;

    #pragma unroll
    for (int oct = 0; oct < 8; ++oct) {
        int childId = childBase + oct;

        int bitX = (oct >> 2) & 1;
        int bitY = (oct >> 1) & 1;
        int bitZ =  oct       & 1;

        center_x[childId] = center_x[parentId] + (bitX ? offset : -offset);
        center_y[childId] = center_y[parentId] + (bitY ? offset : -offset);
        center_z[childId] = center_z[parentId] + (bitZ ? offset : -offset);

        size[childId] = childSize;
        size_sqr[childId] = childSize * childSize;

        parent[childId] = parentId;
        first_child[childId] = -1;

        node_depth[childId] = node_depth[parentId] + 1;
        mass[childId]=0;

        com_x[childId]=0;
        com_y[childId]=0;
        com_z[childId]=0;

        if (childCount[oct] == 0) {
            particle_start[childId] = -1;
            particle_count[childId] = 0;
            continue;
        }

        // If the child is small enough, make it a leaf.
        // Otherwise enqueue it for the next BFS level.
        if (childCount[oct] <= leaf_capacity || depth + 1 >= MORTON_LEVELS - 1) {
            particle_start[childId] = childStart[oct];
            particle_count[childId] = childCount[oct];
        } else {
            particle_start[childId] = -1;
            particle_count[childId] = 0;

            int pos = atomicAdd(d_next_level_count, 1);
            queueFirstNext[pos]   = childStart[oct];
            queueCountNext[pos]   = childCount[oct];
            queueNodeIdxNext[pos] = childId;
        }
    }
}

// Host function: iteratively launch buildOctreeLevelKernel for each tree level.
void BarnesHutForceCalculator::buildOctree(
    const uint32_t* d_sorted_codes,
    int num_particles,
    TreeArrays& tree,
    int* d_next_free_node
) {
    if (num_particles == 0) return;

    // Allocate device queues (size = max possible nodes = num_particles)
    int maxNodes = num_particles;
    int *d_qFirstCur, *d_qCountCur, *d_qIdxCur;
    int *d_qFirstNext, *d_qCountNext, *d_qIdxNext;
    CUDA_CHECK(cudaMalloc(&d_qFirstCur,  maxNodes * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_qCountCur,  maxNodes * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_qIdxCur,    maxNodes * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_qFirstNext, maxNodes * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_qCountNext, maxNodes * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_qIdxNext,   maxNodes * sizeof(int)));

    // Device counter for the next-level queue size
    int *d_next_level_count;
    CUDA_CHECK(cudaMalloc(&d_next_level_count, sizeof(int)));

    // Initialize the queue with the root node [0, num_particles) -> node index 0
    int h_zero = 0, h_np = num_particles;
    CUDA_CHECK(cudaMemcpy(d_qFirstCur, &h_zero, sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_qCountCur, &h_np,   sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_qIdxCur,   &h_zero, sizeof(int), cudaMemcpyHostToDevice));
    int numCurrent = 1;

    // Launch kernels for each depth until no more nodes to split or max depth reached
    for (int depth = 0; depth < constant::MAX_DEPTH && numCurrent > 0; ++depth) {
        // Reset next-level counter to 0
        CUDA_CHECK(cudaMemset(d_next_level_count, 0, sizeof(int)));

        // Launch kernel: one thread per current node
        int blocks = (numCurrent + constant::BLOCK_SIZE - 1) / constant::BLOCK_SIZE;
        buildOctreeLevelKernel<<<blocks, constant::BLOCK_SIZE>>>(
            d_sorted_codes,
            d_qFirstCur, d_qCountCur, d_qIdxCur, numCurrent,
            depth,
            tree.center_x, tree.center_y, tree.center_z,
            tree.size,      tree.size_sqr,
            tree.first_child, tree.parent,
            tree.particle_start, tree.particle_count,
            d_next_free_node,
            d_qFirstNext, d_qCountNext, d_qIdxNext,
            d_next_level_count,
            constant::LEAF_CAPACITY,
            tree.node_depth , tree.mass,
            tree.com_x, tree.com_y, tree.com_z
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        // Read how many children were added
        CUDA_CHECK(cudaMemcpy(&numCurrent, d_next_level_count, sizeof(int), cudaMemcpyDeviceToHost));

        // Swap current and next queues
        std::swap(d_qFirstCur, d_qFirstNext);
        std::swap(d_qCountCur, d_qCountNext);
        std::swap(d_qIdxCur,   d_qIdxNext);
    }

    // Cleanup temporary queues
    CUDA_CHECK(cudaFree(d_qFirstCur));
    CUDA_CHECK(cudaFree(d_qCountCur));
    CUDA_CHECK(cudaFree(d_qIdxCur));
    CUDA_CHECK(cudaFree(d_qFirstNext));
    CUDA_CHECK(cudaFree(d_qCountNext));
    CUDA_CHECK(cudaFree(d_qIdxNext));
    CUDA_CHECK(cudaFree(d_next_level_count));
}

__global__
void computeLeafMassCOMKernel(
    ParticleArrays particles, TreeArrays tree, int num_nodes
)
{
    int node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node >= num_nodes) return;

    // If this is not a leaf, skip
    if (tree.first_child[node] != -1) return;

    int start = tree.particle_start[node];
    int count = tree.particle_count[node];
    if (count <= 0) { tree.mass[node] = 0.0f; return; }

    double m_total = 0.0, com_x = 0.0, com_y = 0.0, com_z = 0.0;
    for(int i = 0; i < count; i++) {
        int p = start + i;
        double m = particles.mass_sorted[p];
        m_total += m;
        com_x += m * particles.pos_x_sorted[p];
        com_y += m * particles.pos_y_sorted[p];
        com_z += m * particles.pos_z_sorted[p];
    }
    tree.mass[node] = static_cast<float>(m_total);
    if (m_total > 0.0) {
        tree.com_x[node] = static_cast<float>(com_x / m_total);
        tree.com_y[node] = static_cast<float>(com_y / m_total);
        tree.com_z[node] = static_cast<float>(com_z / m_total);
    } else {
        tree.com_x[node] = tree.com_y[node] = tree.com_z[node] = 0.0f;
    }
    
}

__global__
void computeInternalMassCOMKernel(
    TreeArrays tree, int num_nodes, int current_depth
) {
    int node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node >= num_nodes) return;

    // Only process nodes at this depth that are internal
    if (tree.node_depth[node] != current_depth) return;
    if (tree.first_child[node] == -1) return;

    double m_total = 0.0, com_x = 0.0, com_y = 0.0, com_z = 0.0;
    int first = tree.first_child[node];

    #pragma unroll
    for (int i = 0; i < 8; i++) {
        int child = first + i;
        if (tree.mass[child] == 0.0f) continue;
        double m = tree.mass[child];
        m_total += m;
        com_x += m * tree.com_x[child];
        com_y += m * tree.com_y[child];
        com_z += m * tree.com_z[child];
    }
    tree.mass[node] = static_cast<float>(m_total);
    if (m_total > 0.0) {
        tree.com_x[node] = static_cast<float>(com_x / m_total);
        tree.com_y[node] = static_cast<float>(com_y / m_total);
        tree.com_z[node] = static_cast<float>(com_z / m_total);
    } else {
        tree.com_x[node] = tree.com_y[node] = tree.com_z[node] = 0.0f;
    }
}


void BarnesHutForceCalculator::computeMassCOM(ParticleArrays& particles) {
    // Get total number of nodes created
    int num_nodes;
    CUDA_CHECK(cudaMemcpy(&num_nodes, nextFreeNode, sizeof(int), cudaMemcpyDeviceToHost));

    int threads = constant::BLOCK_SIZE;
    int blocks = (num_nodes + threads - 1) / threads;

    // First, compute mass and COM for leaf nodes
    computeLeafMassCOMKernel<<<blocks, threads>>>(particles, tree, num_nodes);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Then compute internal nodes bottom-up
    for(int depth = constant::MAX_DEPTH - 1; depth >= 0; --depth) {
        computeInternalMassCOMKernel<<<blocks, threads>>>(tree, num_nodes, depth);
        CUDA_CHECK(cudaDeviceSynchronize());
    }
}


void BarnesHutForceCalculator::buildTree(ParticleArrays& particles) {
    // compute the bounding box of the particles
    computeBoundingBox(particles);

    // Intialize the root node of the tree with the bounding box and all particles
    initializeRootKernel<<<1, 1>>>(d_bounding_cube, tree, nextFreeNode);

    // Compute Morton Encodings for each particle based on their positions and the bounding box
    int threadsPerBlock = 256;
    int blocksPerGrid = (particles.num_particles + constant::BLOCK_SIZE - 1) / constant::BLOCK_SIZE;
    computeMortonCodesKernel<<<blocksPerGrid, threadsPerBlock>>>(particles.pos_x, particles.pos_y, particles.pos_z, d_bounding_cube, particles.morton_codes, particles.particle_indices, particles.num_particles);

    // Radix Sort
    radixSort(particles);

    // sort the particle data based on the sorted indices
    sortParticleDataKernel<<<blocksPerGrid, threadsPerBlock>>>(
        particles.particle_indices_sorted,
        particles.pos_x, particles.pos_y, particles.pos_z, particles.mass,
        particles.pos_x_sorted, particles.pos_y_sorted, particles.pos_z_sorted, particles.mass_sorted,
        particles.num_particles
    );

    // Build the octree using BFS traversal
    buildOctree(particles.morton_codes_sorted, particles.num_particles, tree, nextFreeNode);

    // After building the tree, compute the center of mass and total mass for each node
    computeMassCOM(particles);
}

__global__
void destroyTreeKernel(
    int* next_free_node
)
{
    if (threadIdx.x == 0 && blockIdx.x == 0)
    {
        *next_free_node = 0;
    }
}

void BarnesHutForceCalculator::destroyTree() {
    destroyTreeKernel<<<1, 1>>>(nextFreeNode);
}

void BarnesHutForceCalculator::computeAccelerations(ParticleArrays& particles) {
    if (particles.num_particles == 0) return;

    cudaEventRecord(startEvent, 0);

    // Build tree
    buildTree(particles);

    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    float build_ms = 0;
    cudaEventElapsedTime(&build_ms, startEvent, stopEvent);
    build_tree_time_ms += build_ms;

    cudaEventRecord(startEvent, 0);

    // Traverse tree to compute forces
    int threads = constant::BLOCK_SIZE;
    int blocks = (particles.num_particles + threads - 1) / threads;
    computeAccelerationKernel<<<blocks, threads>>>(
        particles, tree, theta_sqr, softening_squared
    );
    // another kernel which now indexes the acceleration arrays properly
    indexAccelerationKernel<<<blocks, threads>>>(
        particles, particles.acc_x_sorted, particles.acc_y_sorted, particles.acc_z_sorted
    );

    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    float traverse_ms = 0;
    cudaEventElapsedTime(&traverse_ms, startEvent, stopEvent);
    traverse_tree_time_ms += traverse_ms;


    cudaEventRecord(startEvent, 0);

    // destroy the tree after computation
    destroyTree();

    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    float destroy_ms = 0;
    cudaEventElapsedTime(&destroy_ms, startEvent, stopEvent);
    destroy_tree_time_ms += destroy_ms;

    total_compute_time_ms += build_ms + traverse_ms + destroy_ms;
}