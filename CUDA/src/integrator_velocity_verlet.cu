#include "../include/integrator_velocity_verlet.hpp"

__global__
void velocityVerletFirstStepKernel(
    double *pos_x, double *pos_y, double *pos_z,
    double *vel_x, double *vel_y, double *vel_z,
    const double *acc_x, const double *acc_y, const double *acc_z,
    double *old_acc_x, double *old_acc_y, double *old_acc_z,
    int num_particles,
    double dt, double half_dt_sqr
)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_particles) return;
    double local_acc_x = acc_x[tid];
    double local_acc_y = acc_y[tid];
    double local_acc_z = acc_z[tid];

    pos_x[tid] += vel_x[tid] * dt + local_acc_x * half_dt_sqr;
    pos_y[tid] += vel_y[tid] * dt + local_acc_y * half_dt_sqr;
    pos_z[tid] += vel_z[tid] * dt + local_acc_z * half_dt_sqr;

    old_acc_x[tid] = local_acc_x;
    old_acc_y[tid] = local_acc_y;
    old_acc_z[tid] = local_acc_z;
    return;
}

__global__
void velocityVerletSecondStepKernel(
    double *vel_x, double *vel_y, double *vel_z,
    const double *acc_x, const double *acc_y, const double *acc_z,
    const double *old_acc_x, const double *old_acc_y, const double *old_acc_z,
    int num_particles,
    double half_dt
)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_particles) return;
    vel_x[tid] += (old_acc_x[tid] + acc_x[tid]) * half_dt;
    vel_y[tid] += (old_acc_y[tid] + acc_y[tid]) * half_dt;
    vel_z[tid] += (old_acc_z[tid] + acc_z[tid]) * half_dt;
    return;
}

void VelocityVerletIntegrator::step(
    ForceCalculator& calculator,
    ParticleArrays& particles,
    double dt
)
{
    cudaEventRecord(startEvent, 0);
    if(is_first_step) {
        calculator.computeAccelerations(particles);
        is_first_step = false;
    }
    double half_dt = 0.5 * dt;
    double half_dt_sqr = 0.5 * dt * dt;

    int num_particles = particles.num_particles;
    int num_blocks = (num_particles + constant::BLOCK_SIZE - 1) / constant::BLOCK_SIZE;

    velocityVerletFirstStepKernel<<<num_blocks, constant::BLOCK_SIZE>>>(
        particles.pos_x, particles.pos_y, particles.pos_z,
        particles.vel_x, particles.vel_y, particles.vel_z,
        particles.acc_x, particles.acc_y, particles.acc_z,
        particles.old_acc_x, particles.old_acc_y, particles.old_acc_z,
        num_particles,
        dt, half_dt_sqr
    );

    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        std::cerr << "CUDA error in velocityVerletFirstStepKernel: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }

    calculator.computeAccelerations(particles);

    velocityVerletSecondStepKernel<<<num_blocks, constant::BLOCK_SIZE>>>(
        particles.vel_x, particles.vel_y, particles.vel_z,
        particles.acc_x, particles.acc_y, particles.acc_z,
        particles.old_acc_x, particles.old_acc_y, particles.old_acc_z,
        num_particles,
        half_dt
    );

    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    float ms = 0;
    cudaEventElapsedTime(&ms, startEvent, stopEvent);
    total_compute_time_ms += ms;

        cudaDeviceSynchronize();

    err = cudaGetLastError();
    if(err != cudaSuccess) {
        std::cerr << "CUDA error in velocityVerletSecondStepKernel: " << cudaGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
}

std::string VelocityVerletIntegrator::name() const {
    return "Velocity Verlet";
}