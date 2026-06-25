#pragma once
#include "common.hpp"
#include <iostream>

struct Vector3D{
    double x, y, z;
    
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

    // norm of the vector
    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    // normalize the vector
    Vector3D normalize() const {
        double n = norm();
        if (n == 0) return Vector3D(0, 0, 0);
        return Vector3D(x/n, y/n, z/n);
    }

    // Vector Methods and operations
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x+other.x, y+other.y, z+other.z);
    }

    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x-other.x, y-other.y, z-other.z);
    }

    Vector3D operator*(double scalar) const {
        return Vector3D(x*scalar, y*scalar, z*scalar);
    }

    friend Vector3D operator*(double scalar, const Vector3D& vec) {
        return Vector3D(vec.x*scalar, vec.y*scalar, vec.z*scalar);
    }

    Vector3D operator/(double scalar) const {
        return Vector3D(x/scalar, y/scalar, z/scalar);
    }

    Vector3D& operator+=(const Vector3D& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // Dot product
    double dot(const Vector3D& other) const {
        return x*other.x + y*other.y + z*other.z;
    }

    // Cross product
    Vector3D cross(const Vector3D& other) const {
        return Vector3D(y*other.z - z*other.y,
                        z*other.x - x*other.z,
                        x*other.y - y*other.x);
    }
};

struct Particle {
    double mass;
    Vector3D position;
    Vector3D velocity;
    Vector3D acceleration;
    std::string name;
    
    // default constructor
    Particle() : mass(1), position(), velocity(), acceleration(), name("") {}

    // parameterized constructors

    // Constructor with mass, position, velocity with optional name
    Particle(double mass, const Vector3D& position, const Vector3D& velocity, const std::string& name = "")
        : mass(mass), position(position), velocity(velocity), acceleration(), name(name) {}

    // Constructor with mass, position, velocity, acceleration with optional name
    Particle(double mass, const Vector3D& position, const Vector3D& velocity, const Vector3D& acceleration, const std::string& name = "")
        : mass(mass), position(position), velocity(velocity), acceleration(acceleration), name(name) {}

    // function to tell the particle's information
    void info() const {
        if(name.size()) std::cout << "Particle: " << name << "\n";
        std::cout << "Mass: " << mass << "\n";
        std::cout << "\n";
        std::cout << "Position:     (" << position.x << ", " << position.y << ", " << position.z << ")\n";
        std::cout << "Velocity:     (" << velocity.x << ", " << velocity.y << ", " << velocity.z << ")\n";
        std::cout << "Acceleration: (" << acceleration.x << ", " << acceleration.y << ", " << acceleration.z << ")\n";
    }
};