#include "../include/integrator_euler.hpp"

__global__
void eulerIntegrationKernel(
    double *pos_x, double *pos_y, double *pos_z,
    double *vel_x, double *vel_y, double *vel_z,
    const double *acc_x, const double *acc_y, const double *acc_z,
    int num_particles, double dt
)
{
    // okay so in here, we don't need to use shared memory because each thread is only updating its own particle's position and velocity, so we can just use the global memory directly
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= num_particles) return;
    pos_x[tid] += vel_x[tid] * dt;
    pos_y[tid] += vel_y[tid] * dt;
    pos_z[tid] += vel_z[tid] * dt;
    vel_x[tid] += acc_x[tid] * dt;
    vel_y[tid] += acc_y[tid] * dt;
    vel_z[tid] += acc_z[tid] * dt;
    return;
}

void EulerIntegrator::step(ForceCalculator& calculator, ParticleArrays& particles, double dt){
    cudaEventRecord(startEvent, 0); // Start timing

    // compute accelerations using the provided force calculator
    calculator.computeAccelerations(particles);

    // calculate the number of blocks needed to cover all particles
    int num_particles = particles.num_particles;
    int num_blocks = (num_particles + constant::BLOCK_SIZE - 1) / constant::BLOCK_SIZE;

    // launch the kernel to perform the Euler integration step
    eulerIntegrationKernel<<<num_blocks, constant::BLOCK_SIZE>>>(
        particles.pos_x, particles.pos_y, particles.pos_z,
        particles.vel_x, particles.vel_y, particles.vel_z,
        particles.acc_x, particles.acc_y, particles.acc_z,
        num_particles, dt
    );

    cudaEventRecord(stopEvent, 0); // Stop timing
    cudaEventSynchronize(stopEvent);
    float ms = 0;
    cudaEventElapsedTime(&ms, startEvent, stopEvent);
    total_compute_time_ms += ms;
    
    // synchronize the device to ensure all threads have completed before returning
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error in EulerIntegrator::step: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
    cudaDeviceSynchronize();
}

std::string EulerIntegrator::name() const {
    return "Euler Integrator";
}