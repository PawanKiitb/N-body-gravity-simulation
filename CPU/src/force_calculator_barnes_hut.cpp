#include "../include/force_calculator_barnes_hut.hpp"
#include <limits>
#include <functional>

void ForceCalculatorBarnesHut::computeAccelerations(std::vector<Particle>& particles) {
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
    for(Particle &p : particles) {
        p.acceleration = Vector3D(0, 0, 0);
        computeAccelerationForParticle(root, &p);
    }

    end = Clock::now();
    elapsed_time = std::chrono::duration<double, std::milli>(end - start).count();
    traversal.calls++;
    traversal.total_time += elapsed_time;
    total.total_time += elapsed_time;
    // after computing the accelerations, we can destroy the tree to free memory

    start = Clock::now();

    destroyTree(root);
    root = nullptr;
    
    end = Clock::now();
    elapsed_time = std::chrono::duration<double, std::milli>(end - start).count();
    destroyingTree.calls++;
    destroyingTree.total_time += elapsed_time;
    total.total_time += elapsed_time;

    total.calls++;
}

void ForceCalculatorBarnesHut::computeAccelerationForParticle(Node_* node, Particle* particle) {
    if(node == nullptr) return;

    // if node is a leaf node and contains particles, we compute the force directly
    if(node->is_leaf) {
        // now its a leaf node: Two choices: its either the smallest sized node which can contain multiple particles, or its a node with a single particle. In both cases, we can compute the force directly
        for(Particle* p : node->particles) {
            if(p == particle) continue; // skip self-interaction
            Vector3D r = p->position - particle->position;
            double dist_sqr_soft = r.dot(r) + softening_squared;
            double inv_dist = 1.0 / std::sqrt(dist_sqr_soft);
            double inv_dist3 = inv_dist * inv_dist * inv_dist;
            particle->acceleration += constant::G * r * (p->mass * inv_dist3);
        }
        return;
    }
    
    // If node is an internal node, we check the ratio and theta
    else {
        int childIndex = getChildIndex(node, particle->position);
        for(int i = 0; i<8; i++) {
            Node_* child = node->children[i];
            if(child == nullptr) continue; // skip empty children
            if(i == childIndex) {
                // recurse into the child node that contains the particle
                computeAccelerationForParticle(child, particle);
            }
            else {
                // barnes hut condition check:
                Vector3D r = child->center_of_mass - particle->position;
                double dist_sqr = r.dot(r);
                if(child->size_sqr < theta_sqr * dist_sqr) {
                    // treat this child node as a single particle
                    double dist_sqr_soft = dist_sqr + softening_squared;
                    double inv_dist = 1.0 / std::sqrt(dist_sqr_soft);
                    double inv_dist3 = inv_dist * inv_dist * inv_dist;
                    particle->acceleration += constant::G * r * (child->mass * inv_dist3);
                }
                else {
                    // recurse into the child node
                    computeAccelerationForParticle(child, particle);
                }
            }
        }
    } 
}

void ForceCalculatorBarnesHut::buildTree(std::vector<Particle>& particles) {
    if(particles.empty()) {
        root = nullptr;
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

    if(root != nullptr) {
        destroyTree(root);
        root = nullptr;
    }

    root = new Node_();
    root->center = Vector3D((min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2);
    root->size = size;
    root->size_sqr = size * size;
    root->is_leaf = true;

    // Now iterate over all particles
    size_t n = particles.size();
    for(size_t i = 0; i<n; i++) {
        insertParticle(root, &particles[i]);
    }
}

void ForceCalculatorBarnesHut::insertParticle(Node_* node, Particle* particle) {
    // first check if the node is an empty leaf node
    if(node->is_leaf && node->particles.empty()) {
        node->particles.push_back(particle);
        node->mass = particle->mass;
        node->center_of_mass = particle->position;
        return;
    }
    else if(node->is_leaf && !node->particles.empty()) {
        // now we have a leaf node with a particle, we need to subdivide it
        // since its a leaf node, its not subdivided yet, so we need to subdivide it
        if(!subdivideNode(node)) {
            // if we cannot subdivide, we will just put the new particle in the same node, and update the mass and center of mass
            node->mass += particle->mass;
            node->center_of_mass = (node->center_of_mass * (node->mass - particle->mass) + particle->position * particle->mass) / node->mass;
            node->particles.push_back(particle);
            return;
        }
        // subdivison is successful that means that the node only had one particle as multiple particles are only possible in epsilon sized octants, so we can safely assume that the node has only one particle in it, which is the existing particle
        Particle* existing_particle = node->particles.front();
        node->is_leaf = false;
        node->particles.clear(); // clear the particles vector since we are now subdividing

        // now we need to insert the existing particle into the appropriate child node
        int existing_index = getChildIndex(node, existing_particle->position);
        insertParticle(node->children[existing_index], existing_particle);
        // now we can insert the new particle into the appropriate child node
        int new_index = getChildIndex(node, particle->position);
        insertParticle(node->children[new_index], particle);
        // now we need to update the mass and center of mass of the current node
        node->mass = existing_particle->mass + particle->mass;
        node->center_of_mass = (existing_particle->position * existing_particle->mass + particle->position * particle->mass) / node->mass;
    }
    else {
        // if the node is not a leaf, we need to insert the particle into the appropriate child node
        int index = getChildIndex(node, particle->position);
        insertParticle(node->children[index], particle);
        // now we need to update the mass and center of mass of the current node
        node->mass += particle->mass;
        node->center_of_mass = (node->center_of_mass * (node->mass - particle->mass) + particle->position * particle->mass) / node->mass;
    }
}

// Function to determine which child node a particle belongs to which octant, given the node and the particle's position
int ForceCalculatorBarnesHut::getChildIndex(Node_* node, const Vector3D& particle_position) {
    int index = 0;
    if(particle_position.x >= node->center.x) index |= 1; // right half
    if(particle_position.y >= node->center.y) index |= 2; // top half
    if(particle_position.z >= node->center.z) index |= 4; // front half
    return index;
}

// Function to subdivide a node into 8 children, given the node
bool ForceCalculatorBarnesHut::subdivideNode(Node_* node) {
    double half_size = node->size / 2.0;
    if(half_size < epsilon) {
        return false;
    }
    for(int i = 0; i<8; i++) {
        node->children[i] = new Node_();
        node->children[i]->size = half_size;
        node->children[i]->size_sqr = half_size * half_size;
        node->children[i]->is_leaf = true;

        // Determine the center of the child node based on the index
        double offset_x = (i & 1) ? half_size / 2.0 : -half_size / 2.0;
        double offset_y = (i & 2) ? half_size / 2.0 : -half_size / 2.0;
        double offset_z = (i & 4) ? half_size / 2.0 : -half_size / 2.0;

        node->children[i]->center = Vector3D(
            node->center.x + offset_x,
            node->center.y + offset_y,
            node->center.z + offset_z
        );
    }
    return true;
}

void ForceCalculatorBarnesHut::destroyTree(Node_* node) {
    if(node == nullptr) return;
    for(int i = 0; i<8; i++) {
        destroyTree(node->children[i]);
    }
    delete node;
}

std::string ForceCalculatorBarnesHut::name() const {
    return "Barnes-Hut";
}

TimingStatistics& ForceCalculatorBarnesHut::getTimeStatistics() const {
    return const_cast<TimingStatistics&>(total);
}