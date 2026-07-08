#include "../include/force_calculator_naive.hpp"


__global__
void computeAccelerationNaiveKernel(
    const double *pos_x, const double *pos_y, const double *pos_z,
    double *acc_x, double *acc_y, double *acc_z,
    const double *mass, int num_particles,
    double softening_squared
)
{
    __shared__ double shared_pos_x[constant::TILE_SIZE];
    __shared__ double shared_pos_y[constant::TILE_SIZE];
    __shared__ double shared_pos_z[constant::TILE_SIZE];
    __shared__ double shared_mass[constant::TILE_SIZE];

    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const bool valid = tid < num_particles;

    const double my_pos_x = valid ? pos_x[tid] : 0.0;
    const double my_pos_y = valid ? pos_y[tid] : 0.0;
    const double my_pos_z = valid ? pos_z[tid] : 0.0;

    double acc_x_local = 0.0, acc_y_local = 0.0, acc_z_local = 0.0;

    for (int tile_start = 0; tile_start < num_particles; tile_start += constant::TILE_SIZE) {
        int j = tile_start + threadIdx.x;
        if (j < num_particles) {
            shared_pos_x[threadIdx.x] = pos_x[j];
            shared_pos_y[threadIdx.x] = pos_y[j];
            shared_pos_z[threadIdx.x] = pos_z[j];
            shared_mass[threadIdx.x] = mass[j];
        }
        __syncthreads();   // now ALL threads in the block reach this, every time

        if (valid) {
            for (j = 0; j < constant::TILE_SIZE && (tile_start + j) < num_particles; ++j) {
                if (tid == tile_start + j) continue;
                double dx = shared_pos_x[j] - my_pos_x;
                double dy = shared_pos_y[j] - my_pos_y;
                double dz = shared_pos_z[j] - my_pos_z;
                double dist_sqr = dx*dx + dy*dy + dz*dz + softening_squared;
                double inv_dist = 1.0 / sqrt(dist_sqr);
                double inv_dist3 = inv_dist * inv_dist * inv_dist;
                double force = constant::G * shared_mass[j] * inv_dist3;
                acc_x_local += dx * force;
                acc_y_local += dy * force;
                acc_z_local += dz * force;
            }
        }
        __syncthreads();   // ALL threads reach this too
    }

    if (valid) {
        acc_x[tid] = acc_x_local;
        acc_y[tid] = acc_y_local;
        acc_z[tid] = acc_z_local;
    }
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