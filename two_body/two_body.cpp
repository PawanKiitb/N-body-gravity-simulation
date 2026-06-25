#include "common.hpp"
#include "particle.hpp"

triples eccentricity(triples r_rel, triples v_rel, double mu) {
    double r = r_rel.norm();
    double v = v_rel.norm();
    triples h_vec = r_rel ^ v_rel; // angular momentum vector
    return (v_rel ^ h_vec) * (1/mu) - r_rel * (1/r); // eccentricity vector
}

int main() {
    particle a = {};
    particle b = {};

    // set some randome values for the particles
    a.position = {1.0, 0.0, 0.0};
    a.velocity = {0.0, 0.5, 0.0};

    b.position = {-1.0, 0.0, 0.0};
    b.velocity = {0.0, -0.5, 0.0};

    a.mass = b.mass = 1.0;
    
    /************ELLPSES ONLY************/
    /****************START***************/
    triples r_rel = a.position - b.position; // relative position
    triples v_rel = a.velocity - b.velocity; // relative velocity
    double r = r_rel.norm();
    double v = v_rel.norm(); 
    double mu = G*(a.mass + b.mass); 

    triples e_vec = eccentricity(r_rel, v_rel, mu); // eccentricity vector
    double e_mag = e_vec.norm();

    if(e_mag >= 1) {
        std::cerr << "The orbit is not an ellipse (eccentricity >= 1)." << std::endl;
        return 1;
    }

    triples h_vec = r_rel ^ v_rel; // angular momentum vector

    double eps = v*v/2 - mu/r; 
    double a_orbit = -mu/(2*eps);

    triples p_hat = e_vec * (1/e_mag);
    triples h_hat = h_vec * (1/h_vec.norm());

    triples q_hat = h_hat ^ p_hat;
    q_hat = q_hat * (1/q_hat.norm());

    double cosE = (1 - r/a_orbit) / e_mag;
    cosE = std::max(-1.0, std::min(1.0, cosE));
    double E0 = acos(cosE);

    double rdot = (r_rel*v_rel)/r;
    if(rdot < 0) {
        E0 = 2*M_PI - E0;
    }

    double M0 = E0 - e_mag*sin(E0);
    double n = sqrt(mu/(a_orbit*a_orbit*a_orbit));

    std::vector<triples> trajectory_a;
    std::vector<triples> trajectory_b;

    triples R0 = (a.position * a.mass + b.position * b.mass) * (1/(a.mass + b.mass));
    triples V_cm = (a.velocity * a.mass + b.velocity * b.mass) * (1/(a.mass + b.mass));


    std::cout << "eccentricity: " << e_mag << std::endl;
    std::cout << "semi-major axis: " << a_orbit << std::endl;

    for(double t = 0.0; t<10.0; t+=0.1) {
        double M = M0 + n*t;
        double E = M;
        for(int i=0; i<10; i++) {
            E -= (E - e_mag*sin(E) - M) / (1 - e_mag*cos(E));
        }
        double x = a_orbit*(cos(E) - e_mag);
        double y = a_orbit*sqrt(1 - e_mag*e_mag)*sin(E);
        triples r_orbit = p_hat * x + q_hat * y;
        triples R_cm = R0 + V_cm * t;

        triples r_a = R_cm - r_orbit * (b.mass/(a.mass + b.mass));
        triples r_b = R_cm + r_orbit * (a.mass/(a.mass + b.mass));

        trajectory_a.push_back(r_a);
        trajectory_b.push_back(r_b);
    }
    std::ofstream file("orbit.txt");

    for(size_t i = 0; i < trajectory_a.size(); i++) {
        file
            << trajectory_a[i].x << " "
            << trajectory_a[i].y << " "
            << trajectory_a[i].z << " "
            << trajectory_b[i].x << " "
            << trajectory_b[i].y << " "
            << trajectory_b[i].z << "\n";
    }

    file.close();
}