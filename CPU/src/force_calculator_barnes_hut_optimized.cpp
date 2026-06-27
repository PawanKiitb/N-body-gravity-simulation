#include "../include/force_calculator_barnes_hut_optimized.hpp"
#include <limits>
#include <functional>
#include <algorithm>

void ForceCalculatorBarnesHutOptimized::computeAccelerations(std::vector<Particle>& particles) {
    particles_ = &particles;
    // firstly build the tree
    auto start = Clock::now();

    buildTree(particles);

    auto end = Clock::now();
    double elapsed_time = std::chrono::duration<double, std::milli>(end - start).count();
    buildingTree.calls++;
    buildingTree.total_time += elapsed_time;
    total.total_time += elapsed_time;

    // so the tree is built, now we need to compute the accelerations for each particle
    start = Clock::now();

    #pragma omp parallel for schedule(static)
    for(size_t i = 0; i<particles.size(); i++) {
        particles[i].acceleration = Vector3D(0, 0, 0);
        computeAccelerationForParticle(rootIndex, i);
    }

    end = Clock::now();
    elapsed_time = std::chrono::duration<double, std::milli>(end - start).count();
    traversal.calls++;
    traversal.total_time += elapsed_time;
    total.total_time += elapsed_time;
    // after computing the accelerations, we can destroy the tree to free memory

    start = Clock::now();

    destroyTree();
    
    end = Clock::now();
    elapsed_time = std::chrono::duration<double, std::milli>(end - start).count();
    destroyingTree.calls++;
    destroyingTree.total_time += elapsed_time;
    total.total_time += elapsed_time;

    total.calls++;
}

void ForceCalculatorBarnesHutOptimized::computeAccelerationForParticle(int rootIndex, int particleIndex) {
    if(rootIndex == -1) return;
    Particle& particle = (*particles_)[particleIndex];
    // Maximum tree depth is very small (typically < 32)
    int stack[64];
    int top = 0;
    stack[top++] = rootIndex;
    while(top > 0) {
        int nodeIndex = stack[--top];
        Node& node = nodePool[nodeIndex];
        // Leaf node
        if(node.firstChild == -1) {
            for(size_t i = 0; i < node.particleIndices.size(); i++) {
                int otherIndex = node.particleIndices[i];
                if(otherIndex == particleIndex) continue;
                const Particle& other = (*particles_)[otherIndex];
                Vector3D r = other.position - particle.position;
                double distSqrSoft = r.dot(r) + softening_squared;
                double invDist = 1.0 / std::sqrt(distSqrSoft);
                double invDist3 = invDist * invDist * invDist;
                particle.acceleration += constant::G * r * (other.mass * invDist3);
            }
            continue;
        }
        // Internal node
        // int containingChild = getChildIndex(nodeIndex, particle.position);
        int containingChild = 0;
        if(particle.position.x >= node.center.x) containingChild |= 1;
        if(particle.position.y >= node.center.y) containingChild |= 2;
        if(particle.position.z >= node.center.z) containingChild |= 4;
        Node* children = &nodePool[node.firstChild];
        for(int i = 0; i < 8; i++) {
            Node& child = children[i];
            if(i == containingChild) {
                stack[top++] = node.firstChild + i;
                continue;
            }
            Vector3D r = child.center_of_mass - particle.position;
            double distSqr = r.dot(r);
            if(child.size_sqr < theta_sqr * distSqr) {
                double distSqrSoft = distSqr + softening_squared;
                double invDist = 1.0 / std::sqrt(distSqrSoft);
                double invDist3 = invDist * invDist * invDist;
                particle.acceleration += constant::G * r * (child.mass * invDist3);
            }
            else {
                stack[top++] = node.firstChild + i;
            }
        }
    }
}
void ForceCalculatorBarnesHutOptimized::buildTree(std::vector<Particle>& particles) {
    if(particles.empty()) {
        rootIndex = -1;
        return;
    }
    // first find the bounding box of all particles
    double min_x, max_x, min_y, max_y, min_z, max_z;
    min_x = min_y = min_z = std::numeric_limits<double>::max();
    max_x = max_y = max_z = std::numeric_limits<double>::lowest();
    for(Particle &p : particles) {
        min_x = std::min(min_x, p.position.x);
        max_x = std::max(max_x, p.position.x);
        min_y = std::min(min_y, p.position.y);
        max_y = std::max(max_y, p.position.y);
        min_z = std::min(min_z, p.position.z);
        max_z = std::max(max_z, p.position.z);
    }

    double size = std::max(max_x - min_x, std::max(max_y - min_y, max_z - min_z)); 
    size += epsilon;

    rootIndex = allocateNode();
    Node& root = nodePool[rootIndex];
    root.center = Vector3D((min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2);
    root.size = size;
    root.size_sqr = size * size;

    for(size_t i = 0; i<particles.size(); i++) {
        insertParticle(rootIndex, i);
    }
}

void ForceCalculatorBarnesHutOptimized::insertParticle(int nodeIndex, int particleIndex) {
    Node& node = nodePool[nodeIndex];
    Particle& particle = (*particles_)[particleIndex];
    node.mass += particle.mass;
    node.center_of_mass = (node.center_of_mass * (node.mass - particle.mass) + particle.position * particle.mass) / node.mass;

    // If the node is a leaf
    if(node.firstChild == -1) {
        if(node.particleIndices.size() < leaf_capacity) {
            // that means this node is a leaf and has space for more particles.
            // So we simply add the particle to this node
            node.particleIndices.push_back(particleIndex);
            return;
        }
        else {
            // that means this node is a leaf but is full, so we need to subdivide it
            if(!subdivideNode(nodeIndex)) {
                // if we cannot subdivide, we simply add the particle to this node
                node.particleIndices.push_back(particleIndex);
                return;
            }
            // now the subdivion was successful, we need to reinsert the existing particles into the children
            std::vector<int> oldParticleIndices;
            oldParticleIndices.swap(node.particleIndices);
            for(size_t i = 0; i<oldParticleIndices.size(); i++) {
                int existingParticleIndex = oldParticleIndices[i];
                const Vector3D& pos = (*particles_)[existingParticleIndex].position;
                int childIndex = 0;
                if(pos.x >= node.center.x) childIndex |= 1;
                if(pos.y >= node.center.y) childIndex |= 2;
                if(pos.z >= node.center.z) childIndex |= 4;
                insertParticle(node.firstChild + childIndex, existingParticleIndex);
            }
            // now we can insert the new particle into the appropriate child
            // int childIndex = getChildIndex(nodeIndex, particle.position);
            int childIndex = 0;
            if(particle.position.x >= node.center.x) childIndex |= 1;
            if(particle.position.y >= node.center.y) childIndex |= 2;
            if(particle.position.z >= node.center.z) childIndex |= 4;
            insertParticle(node.firstChild + childIndex, particleIndex);
            // node.particleIndices.clear(); // clear the particle indices since this node is no longer a leaf
            return;
        }
    }
    else {
        // that means that this node is an internal node, so we need to insert the particle into the appropriate child
        // int childIndex = getChildIndex(nodeIndex, particle.position);
        int childIndex = 0;
        if(particle.position.x >= node.center.x) childIndex |= 1;
        if(particle.position.y >= node.center.y) childIndex |= 2;
        if(particle.position.z >= node.center.z) childIndex |= 4;
        insertParticle(node.firstChild + childIndex, particleIndex);
        return;
    }
}

// Allocate a new node
int ForceCalculatorBarnesHutOptimized::allocateNode() {
    assert(nextFreeNode < nodePool.size());
    Node* node = &nodePool[nextFreeNode++];
    *node = Node(); // Reset the node
    return nextFreeNode - 1;
}

// Function to determine which child node a particle belongs to which octant, given the node and the particle's position
// int ForceCalculatorBarnesHutOptimized::getChildIndex(int nodeIndex, const Vector3D& particle_position) {
//     int index = 0;
//     Node& node = nodePool[nodeIndex];
//     if(particle_position.x >= node.center.x) index |= 1;
//     if(particle_position.y >= node.center.y) index |= 2;
//     if(particle_position.z >= node.center.z) index |= 4;
//     return index;
// }

// Function to subdivide a node into 8 children, given the node
bool ForceCalculatorBarnesHutOptimized::subdivideNode(int nodeIndex) {
    Node& node = nodePool[nodeIndex];
    double half_size = node.size / 2.0;
    if(half_size < epsilon) {
        return false;
    }
    node.firstChild = nextFreeNode;
    // node.is_leaf = false;
    for(int i = 0; i<8; i++) {
        Node& child = nodePool[allocateNode()];
        child.size = half_size;
        child.size_sqr = half_size * half_size;

        double offset_x = (i & 1) ? half_size / 2.0 : -half_size / 2.0;
        double offset_y = (i & 2) ? half_size / 2.0 : -half_size / 2.0;
        double offset_z = (i & 4) ? half_size / 2.0 : -half_size / 2.0;

        child.center = Vector3D(
            node.center.x + offset_x,
            node.center.y + offset_y,
            node.center.z + offset_z
        );
        // child.is_leaf = true;
    }
    return true;
}

void ForceCalculatorBarnesHutOptimized::destroyTree() {
    nextFreeNode = 0;
    rootIndex = -1;
}

std::string ForceCalculatorBarnesHutOptimized::name() const {
    return "Barnes-Hut";
}

TimingStatistics& ForceCalculatorBarnesHutOptimized::getTimeStatistics() const {
    return const_cast<TimingStatistics&>(total);
}