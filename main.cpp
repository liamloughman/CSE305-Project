#include <iostream>
#include <vector>
#include <cmath>

const double G = 6.67430e-11;

struct Body {
    double mass;
    double x, y;
    double vx, vy;
};

void updateForcesAndPositions(std::vector<Body>& bodies, double dt) {
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double distSq = dx * dx + dy * dy;
            double dist = std::sqrt(distSq);
            double force = G * bodies[i].mass * bodies[j].mass / distSq;
            double fx_i = force * dx / dist;
            double fy_i = force * dy / dist;
            fx[i] += fx_i;
            fy[i] += fy_i;
            fx[j] -= fx_i;
            fy[j] -= fy_i;
        }
    }

    for (int i = 0; i < n; ++i) {
        bodies[i].vx += fx[i] / bodies[i].mass * dt;
        bodies[i].vy += fy[i] / bodies[i].mass * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
}

int main() {
    
    std::vector<Body> bodies = {
        {1e10, 0, 0, 0, 0}, // mass, x, y, vx, vy
        {5e10, 1e5, 0, 0, -10}
    };

    double dt = 1;  // time step 
    int steps = 3;  // total number of steps

    for (int step = 0; step < steps; ++step) {
        updateForcesAndPositions(bodies, dt);
        std::cout << "Step " << step << ":\n";
        for (const Body& body : bodies) {
            std::cout << "Body at (" << body.x << ", " << body.y << ")\n";
        }
    }
    return 0;
}