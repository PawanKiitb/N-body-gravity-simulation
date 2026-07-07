#include "../src/test_common.cpp"

#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

int main() {
    std::ifstream file("/home/pawan/Sem4/SOC_2026/CPU/tests/particles_100.txt");
    if(!file.is_open()) {
        std::cerr << "Error opening file particles_100.txt\n";
        return 1;
    }
    // file format:
    // N
    // m x y z vx vy vz
    std::vector<double> h_pos_x;
    std::vector<double> h_pos_y;
    std::vector<double> h_pos_z;
    std::vector<double> h_vel_x;
    std::vector<double> h_vel_y;
    std::vector<double> h_vel_z;
    std::vector<double> h_mass;
    int N;
    file >> N;
    h_pos_x.resize(N);
    h_pos_y.resize(N);
    h_pos_z.resize(N);
    h_vel_x.resize(N);
    h_vel_y.resize(N);
    h_vel_z.resize(N);
    h_mass.resize(N);
    for(int i = 0; i < N; i++) {
        file >> h_mass[i] >> h_pos_x[i] >> h_pos_y[i] >> h_pos_z[i]
             >> h_vel_x[i] >> h_vel_y[i] >> h_vel_z[i];
    }
    file.close();

    ParticleArrays bhParticles;
    allocateParticles(bhParticles, N);
    copyParticlesToDevice(
        bhParticles,
        h_pos_x.data(),
        h_pos_y.data(),
        h_pos_z.data(),
        h_vel_x.data(),
        h_vel_y.data(),
        h_vel_z.data(),
        h_mass.data()
    );
    BarnesHutForceCalculator bh(0.0, 0.5);
    bh.computeAccelerations(bhParticles);
}