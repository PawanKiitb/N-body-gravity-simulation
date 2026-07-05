#include "../include/integrator_leapfrog.hpp"

__global__
void LeafFrogFirstStepKernel(
    double *pos_x, double *pos_y, double *pos_z,
    double *vel_x, double *vel_y, double *vel_z,
    const double *acc_x, const double *acc_y, const double *acc_z,
    int num_particles,
    double half_dt, double dt
)
{
    // declare the thread id for the current thread
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= num_particles) return;
    double local_vel_x = vel_x[tid] + acc_x[tid] * half_dt;
    double local_vel_y = vel_y[tid] + acc_y[tid] * half_dt;
    double local_vel_z = vel_z[tid] + acc_z[tid] * half_dt;

    pos_x[tid] += local_vel_x * dt;
    pos_y[tid] += local_vel_y * dt;
    pos_z[tid] += local_vel_z * dt;

    vel_x[tid] = local_vel_x;
    vel_y[tid] = local_vel_y;
    vel_z[tid] = local_vel_z;

    return;
}

__global__
void LeapfrogSecondStepKernel(
    double *vel_x, double *vel_y, double *vel_z,
    const double *acc_x, const double *acc_y, const double *acc_z,
    int num_particles,
    double half_dt
)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_particles) return;
    vel_x[tid] += acc_x[tid] * half_dt;
    vel_y[tid] += acc_y[tid] * half_dt;
    vel_z[tid] += acc_z[tid] * half_dt;
    return;
}

void LeapfrogIntegrator::step(
    ForceCalculator& calculator,
    ParticleArrays& particles,
    double dt
)
{
    cudaEventRecord(startEvent, 0); // Start timing

    // get the number of particles
    int num_particles = particles.num_particles;
    double half_dt = 0.5 * dt;

    if(is_first_step) {
        calculator.computeAccelerations(particles);
        is_first_step = false;
    }
    // calculate the number of blocks needed to cover all particles
    int num_blocks = (num_particles + constant::BLOCK_SIZE - 1) / constant::BLOCK_SIZE;
    
    LeapfrogFirstStepKernel<<<num_blocks, constant::BLOCK_SIZE>>>(
        particles.pos_x, particles.pos_y, particles.pos_z,
        particles.vel_x, particles.vel_y, particles.vel_z,
        particles.acc_x, particles.acc_y, particles.acc_z,
        num_particles,
        half_dt, dt
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error in LeapfrogIntegrator::step: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
    cudaDeviceSynchronize();

    calculator.computeAccelerations(particles);

    LeapfrogSecondStepKernel<<<num_blocks, constant::BLOCK_SIZE>>>(
        particles.vel_x, particles.vel_y, particles.vel_z,
        particles.acc_x, particles.acc_y, particles.acc_z,
        num_particles,
        half_dt
    );

    cudaEventRecord(stopEvent, 0); // Stop timing
    cudaEventSynchronize(stopEvent);
    float ms = 0;
    cudaEventElapsedTime(&ms, startEvent, stopEvent);
    total_compute_time_ms += ms;

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error in LeapfrogIntegrator::step: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
    cudaDeviceSynchronize();

}

std::string LeapfrogIntegrator::name() const {
    return "Leapfrog Integrator";
}