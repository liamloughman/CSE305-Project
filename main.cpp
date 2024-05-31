#include <iostream>
#include <vector>
#include <cmath>
#include <Magick++.h>
#include <thread>
#include <mutex>

const double G = 6.67430e-11;

struct Body {
    double mass;
    double x, y;
    double vx, vy;
};


std::mutex mutex;

void updateForces(std::vector<Body>& bodies, std::vector<double>& fx, std::vector<double>& fy) {
    int n = bodies.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
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
                fy[j] -= fy_i;}
        }
    }

}

void Seq(std::vector<Body>& bodies, double dt) {
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);
    updateForces(bodies, fx, fy);

    for (int i = 0; i < n; ++i) {
        bodies[i].vx += fx[i] / bodies[i].mass * dt;
        bodies[i].vy += fy[i] / bodies[i].mass * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
}

void updateVelocities(std::vector<Body>& bodies, double dt, std::vector<double>& fx, std::vector<double>& fy, int start, int end) {
    for (int i = start; i < end; ++i) {
        bodies[i].vx += fx[i] / bodies[i].mass * dt;
        bodies[i].vy += fy[i] / bodies[i].mass * dt;
    }
}

void updatePositions(std::vector<Body>& bodies, double dt, int start, int end) {
    for (int i = start; i < end; ++i) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
}

void Step_Parallel(std::vector<Body>& bodies, double dt) {
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);
    updateForces(bodies, fx, fy);

    int numThreads = std::thread::hardware_concurrency();
    int blockSize = n / numThreads;
    int remainder = n % numThreads;

    std::vector<std::thread> threads1(numThreads);

    for (int i = 0; i < numThreads; ++i) {
        int start = i * blockSize + std::min(i, remainder);
        int end = start + blockSize + (i < remainder ? 1 : 0);
        threads1[i] = std::thread(updateVelocities, std::ref(bodies), dt, std::ref(fx), std::ref(fy), start, end);
    }

    for (std::thread &thread : threads1) {
        thread.join();
    }

    std::vector<std::thread> threads2(numThreads);

    for (int i = 0; i < numThreads; ++i) {
        int start = i * blockSize + std::min(i, remainder);
        int end = start + blockSize + (i < remainder ? 1 : 0);
        threads2[i] = std::thread(updatePositions, std::ref(bodies), dt, start, end);
    }

    for (std::thread &thread : threads2) {
        thread.join();
    }
}

void updateForcesParallel(std::vector<Body>& bodies, std::vector<double>& fx, std::vector<double>& fy, int start, int end) {
    int n = bodies.size();
    for (int i = start; i < end; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double distSq = dx * dx + dy * dy;
                double dist = std::sqrt(distSq);
                double force = G * bodies[i].mass * bodies[j].mass / distSq;
                double fx_i = force * dx / dist;
                double fy_i = force * dy / dist;

                mutex.lock();
                fx[i] += fx_i;
                fy[i] += fy_i;
                fx[j] -= fx_i;
                fy[j] -= fy_i;
                mutex.unlock();
            }
        }
    }
}

void ForcesParallel(std::vector<Body>& bodies, double dt) {
    int numThreads = std::thread::hardware_concurrency();
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);

    std::vector<std::thread> threads;

    for (int i = 0; i < numThreads; ++i) {
        int start = i * n / numThreads;
        int end = (i + 1) * n / numThreads;
        threads.push_back(std::thread(updateForcesParallel, std::ref(bodies), std::ref(fx), std::ref(fy), start, end));
    }
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();

    for (int i = 0; i < numThreads; ++i) {
        int start = i * n / numThreads;
        int end = (i + 1) * n / numThreads;
        threads.push_back(std::thread(updateVelocities, std::ref(bodies), dt, std::ref(fx), std::ref(fy), start, end));
    }
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();

    for (int i = 0; i < numThreads; ++i) {
        int start = i * n / numThreads;
        int end = (i + 1) * n / numThreads;
        threads.push_back(std::thread(updatePositions, std::ref(bodies), dt, start, end));
    }
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();
}

Magick::Image drawFrame(const std::vector<Body>& bodies) {
    Magick::Image image(Magick::Geometry(800, 600), "black");
    image.type(Magick::TrueColorType);
    image.strokeColor("white");
    image.fillColor("white");

    for (const Body& body : bodies) {
        double x = body.x / 1e9 + 400;
        double y = body.y / 1e9 + 300;
        image.draw(Magick::DrawableCircle(x, y, x + 2, y + 2));
    }
    return image;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <1 for Seq, 2 for Step_Parallel>" << std::endl;
        return 1;
    }

    int algo = std::atoi(argv[1]);
    if (algo != 1 && algo != 2 && algo != 3) {
        std::cerr << "Invalid option. Use 1 for Seq, 2 for Step_Parallel." << std::endl;
        return 1;
    }
    
    std::vector<Body> bodies = {
        {1e25, 0, 0, 0, 0}, // mass, x, y, vx, vy
        {5e2, 5e10, 0, -150, 80}
    };

    double dt = 1e7;  // time step in seconds
    int steps = 400;  // total number of steps

    std::vector<Magick::Image> frames;
    for (int step = 0; step < steps; ++step) {
        if (algo == 1) {
            Seq(bodies, dt);
        } else if (algo == 2){
            Step_Parallel(bodies, dt);
        }
        else if (algo == 3) {
            ForcesParallel(bodies, dt);
        }
        frames.push_back(drawFrame(bodies));
    }

    for (auto& frame : frames) {
        frame.animationDelay(2);
    }
    Magick::writeImages(frames.begin(), frames.end(), "simulation.gif");

    return 0;
}