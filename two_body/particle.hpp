#include <cmath>

struct triples {
    double x, y, z;

    // overload the + operator for position
    triples operator+(const triples& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    // overload the - operator for position
    triples operator-(const triples& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    // overload the * operator for position
    triples operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    // overload the dot product operator for position
    double operator*(const triples& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // overload the cross product operator for position
    triples operator^(const triples& other) const {
        return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
    }

    // get the functionality for norm
    double norm() const {
        return sqrt(x*x + y*y + z*z);
    }
};

struct particle {
    triples position;
    triples velocity;
    double mass;
};