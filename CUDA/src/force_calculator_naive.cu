#include "../include/force_calculator_naive.hpp"

__global__ 
void computeAccelerationNaiveKernel(
    const double *pos_x, const double *pos_y, const double *pos_z,
    double *acc_x, double *acc_y, double *acc_z,
    const double *mass, int num_particles,
    double softening_squared
) 
{
    // firstly, we declare shared memory arrays to hold the positions and masses of particles in the current tile
    __shared__ double shared_pos_x[constant::TILE_SIZE];
    __shared__ double shared_pos_y[constant::TILE_SIZE];
    __shared__ double shared_pos_z[constant::TILE_SIZE];
    __shared__ double shared_mass[constant::TILE_SIZE];
    
    // declare the thread id for the current thread
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // if the thread id is greater than or equal to the number of particles, we return from the kernel
    if(tid >= num_particles) return;
    
    // load the position of the current particle into local variables
    // Each thread will call its own particle's position and mass from global memory
    // we store the position of the current particle in local variables(registers) for faster access
    const double my_pos_x = pos_x[tid];
    const double my_pos_y = pos_y[tid];
    const double my_pos_z = pos_z[tid];

    double acc_x_local = 0.0;
    double acc_y_local = 0.0;
    double acc_z_local = 0.0;

    for(int tile_start = 0; tile_start < num_particles; tile_start += constant::TILE_SIZE) {
        // is the current thread's index within the tile?
        int j = tile_start + threadIdx.x;

        // if it is, we load the position and mass of the particle into shared memory
        if(j < num_particles) {
            shared_pos_x[threadIdx.x] = pos_x[j];
            shared_pos_y[threadIdx.x] = pos_y[j];
            shared_pos_z[threadIdx.x] = pos_z[j];
            shared_mass[threadIdx.x] = mass[j];
        }
        // let all threads come here and wait for the shared memory to be populated before proceeding
        __syncthreads();

        // now we are done with loading the shared memory, we can compute the acceleration for the current particle from the particles in the current tile
        for(j = 0; j < constant::TILE_SIZE && (tile_start + j) < num_particles; ++j) {
            if(tid == tile_start + j) continue; // skip self-interaction
            double dx = shared_pos_x[j] - my_pos_x;
            double dy = shared_pos_y[j] - my_pos_y;
            double dz = shared_pos_z[j] - my_pos_z;
            double dist_sqr = dx * dx + dy * dy + dz * dz + softening_squared;
            double inv_dist = 1.0 / sqrt(dist_sqr);
            double inv_dist3 = inv_dist * inv_dist * inv_dist;
            double force = constant::G * shared_mass[j] * inv_dist3;
            acc_x_local += dx * force;
            acc_y_local += dy * force;
            acc_z_local += dz * force;
        }
        // let all threads come here and wait for the shared memory to be populated before proceeding to the next tile
        __syncthreads();
    }

    // store the computed acceleration in global memory
    acc_x[tid] = acc_x_local;
    acc_y[tid] = acc_y_local;
    acc_z[tid] = acc_z_local;
    return;
}

void ForceCalculatorNaive::computeAccelerations(ParticleArrays& particles) {
    // get the number of particles
    int num_particles = particles.num_particles;

    // calculate the number of blocks needed to cover all particles
    int num_blocks = (num_particles + constant::BLOCK_SIZE - 1) / constant::BLOCK_SIZE;

    cudaEventRecord(startEvent, 0); // Start timing

    // launch the kernel to compute accelerations
    computeAccelerationNaiveKernel<<<num_blocks, constant::BLOCK_SIZE>>>(
        particles.pos_x, particles.pos_y, particles.pos_z,
        particles.acc_x, particles.acc_y, particles.acc_z,
        particles.mass, num_particles,
        softening_squared
    );

    cudaEventRecord(stopEvent, 0); // Stop timing
    cudaEventSynchronize(stopEvent);
    float ms = 0;
    cudaEventElapsedTime(&ms, startEvent, stopEvent);
    total_compute_time_ms += ms;

    // synchronize the device to ensure all threads have completed before returning
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error in computeAccelerations: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }

}

std::string ForceCalculatorNaive::name() const {
    return "Naive-CUDA";
}