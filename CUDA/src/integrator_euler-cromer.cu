#include "../include/integrator_euler-cromer.hpp"

__global__
void eulerCromerIntegrationKernel(
    double *pos_x, double *pos_y, double *pos_z,
    double *vel_x, double *vel_y, double *vel_z,
    const double *acc_x, const double *acc_y, const double *acc_z,
    int num_particles, double dt
)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= num_particles) return;
    double new_vel_x = vel_x[tid] + acc_x[tid] * dt;
    double new_vel_y = vel_y[tid] + acc_y[tid] * dt;
    double new_vel_z = vel_z[tid] + acc_z[tid] * dt;
    pos_x[tid] += new_vel_x * dt;
    pos_y[tid] += new_vel_y * dt;
    pos_z[tid] += new_vel_z * dt;
    vel_x[tid] = new_vel_x;
    vel_y[tid] = new_vel_y;
    vel_z[tid] = new_vel_z;
    
    return;
}

void EulerCromerIntegrator::step(ForceCalculator& calculator, ParticleArrays& particles, double dt){
    cudaEventRecord(startEvent, 0); // Start timing

    // compute accelerations using the provided force calculator
    calculator.computeAccelerations(particles);

    // calculate the number of blocks needed to cover all particles
    int num_particles = particles.num_particles;
    int num_blocks = (num_particles + constant::BLOCK_SIZE - 1) / constant::BLOCK_SIZE;

    // launch the kernel to perform the Euler integration step
    eulerCromerIntegrationKernel<<<num_blocks, constant::BLOCK_SIZE>>>(
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
        std::cerr << "CUDA error in EulerCromerIntegrator::step: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
    cudaDeviceSynchronize();
}

std::string EulerCromerIntegrator::name() const {
    return "Euler-Cromer Integrator";
}