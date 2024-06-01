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

void compute_forces(std::vector<Body>& bodies, std::vector<double>& fx, std::vector<double>& fy) {
    int n = bodies.size();
    for (int i = 0; i < n; ++i) {
        fx[i] = 0;
        fy[i] = 0;
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

void sequential_simulation(std::vector<Body>& bodies, double dt) {
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);
    compute_forces(bodies, fx, fy);

    for (int i = 0; i < n; ++i) {
        bodies[i].vx += fx[i] / bodies[i].mass * dt;
        bodies[i].vy += fy[i] / bodies[i].mass * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
}

void compute_velocities(std::vector<Body>& bodies, double dt, std::vector<double>& fx, std::vector<double>& fy, int start, int end) {
    for (int i = start; i < end; ++i) {
        bodies[i].vx += fx[i] / bodies[i].mass * dt;
        bodies[i].vy += fy[i] / bodies[i].mass * dt;
    }
}

void compute_positions(std::vector<Body>& bodies, double dt, int start, int end) {
    for (int i = start; i < end; ++i) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
}

void parallel_step_simulation(std::vector<Body>& bodies, double dt) {
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);
    compute_forces(bodies, fx, fy);

    int numThreads = std::thread::hardware_concurrency();
    int blockSize = n / numThreads;
    int remainder = n % numThreads;

    std::vector<std::thread> threads1(numThreads);

    for (int i = 0; i < numThreads; ++i) {
        int start = i * blockSize + std::min(i, remainder);
        int end = start + blockSize + (i < remainder ? 1 : 0);
        threads1[i] = std::thread(compute_velocities, std::ref(bodies), dt, std::ref(fx), std::ref(fy), start, end);
    }

    for (std::thread &thread : threads1) {
        thread.join();
    }

    std::vector<std::thread> threads2(numThreads);

    for (int i = 0; i < numThreads; ++i) {
        int start = i * blockSize + std::min(i, remainder);
        int end = start + blockSize + (i < remainder ? 1 : 0);
        threads2[i] = std::thread(compute_positions, std::ref(bodies), dt, start, end);
    }

    for (std::thread &thread : threads2) {
        thread.join();
    }
}

void parallel_distinc_aux(std::vector<Body>& bodies, std::vector<double>& fx, std::vector<double>& fy, int start, int end) {
    int n = bodies.size();
    for (int i = start; i < end; ++i) {
        fx[i] = 0;
        fy[i] = 0;
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

void parallel_distinc_simulation(std::vector<Body>& bodies, double dt) {
    int numThreads = std::thread::hardware_concurrency();
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);

    std::vector<std::thread> threads;

    for (int i = 0; i < numThreads; ++i) {
        int start = i * n / numThreads;
        int end = (i + 1) * n / numThreads;
        threads.push_back(std::thread(parallel_distinc_aux, std::ref(bodies), std::ref(fx), std::ref(fy), start, end));
    }
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();

    for (int i = 0; i < numThreads; ++i) {
        int start = i * n / numThreads;
        int end = (i + 1) * n / numThreads;
        threads.push_back(std::thread(compute_velocities, std::ref(bodies), dt, std::ref(fx), std::ref(fy), start, end));
    }
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();

    for (int i = 0; i < numThreads; ++i) {
        int start = i * n / numThreads;
        int end = (i + 1) * n / numThreads;
        threads.push_back(std::thread(compute_positions, std::ref(bodies), dt, start, end));
    }
    for (auto& thread : threads) {
        thread.join();
    }
    threads.clear();
}

void parallel_combined_aux(std::vector<Body>& bodies, std::vector<Body>& curr, std::vector<Body>& next, size_t start, size_t end, double dt) {
    size_t n = bodies.size();
    for (size_t i = start; i < end; i++) {
        double net_fx = 0;
        double net_fy = 0;
        for (size_t j = 0; j < n; j++) {
            if (i != j) {
                double dx = curr[j].x - curr[i].x;
                double dy = curr[j].y - curr[i].y;
                double dist_sq = dx * dx + dy * dy;
                double dist = std::sqrt(dist_sq);
                double force_mag = G * bodies[i].mass * bodies[j].mass / dist_sq;
                double fx = force_mag * dx / dist;
                double fy = force_mag * dy / dist;
                net_fx += fx;
                net_fy += fy;
            }
        }
        bodies[i].vx += net_fx * (dt / bodies[i].mass);
        bodies[i].vy += net_fy * (dt / bodies[i].mass);
        next[i].x = curr[i].x + bodies[i].vx * dt;
        next[i].y = curr[i].y + bodies[i].vy * dt;
    }
}

void parallel_combined_simulation(std::vector<Body>& bodies, double dt) {
    size_t numThreads = std::thread::hardware_concurrency();
    size_t n = bodies.size();
    std::vector<std::thread> threads(numThreads);
    std::vector<Body> tmp_bodies(n), curr = bodies, next = tmp_bodies;

    size_t chunk_size = (n + numThreads - 1) / numThreads;
    for (size_t i = 0; i < numThreads; ++i) {
        size_t start = i * n / numThreads;
        size_t end = (i + 1) * n / numThreads;
        threads[i] = std::thread(parallel_combined_aux, std::ref(bodies), std::ref(curr), std::ref(next), start, end, dt);
    }
    for (auto& thread : threads) {
        thread.join();
    }
    for (size_t i = 0; i < n; i++) {
        std::swap(curr[i].x, next[i].x);
        std::swap(curr[i].y, next[i].y);
    }
    
    for (size_t i = 0; i < n; i++) {
        std::swap(bodies[i].x, curr[i].x);
        std::swap(bodies[i].y, curr[i].y);
    }
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

std::vector<Body> generate_random_bodies(int N) {
    std::vector<Body> bodies;
    bodies.reserve(N);
    std::srand(std::time(nullptr));

    for (int i = 0; i < N; ++i) {
        double mass = 1e2 + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (1e24 - 1e2)));
        double x = -5e10 + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (5e10 - (-5e10))));
        double y = -5e10 + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (5e10 - (-5e10))));
        double vx = -100 + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (100 - (-100))));
        double vy = -100 + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (100 - (-100))));
        bodies.push_back({mass, x, y, vx, vy});
    }
    return bodies;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <1 for sequential_simulation, 2 for parallel_step_simulation, 3 for parallel_distinc_simulation, 4 for parallel_combined_simulation>" << std::endl;
        return 1;
    }

    int algo = std::atoi(argv[1]);
    if (algo != 1 && algo != 2 && algo != 3 && algo != 4) {
        std::cerr << "Invalid option. Use 1 for sequential_simulation, 2 for parallel_step_simulation, 3 for parallel_distinc_simulation, 4 for parallel_combined_simulation." << std::endl;
        return 1;
    }
    
    /*std::vector<Body> bodies = {
        {1e24, 0, 0, 0, 0}, // mass, x, y, vx, vy
        {1e2, 5e10, 1e10, -70, 60},
        {1e2, -5e10, 5e10, 10, -20}
    };*/

    double dt = 1e7;  // time step in seconds
    int steps = 200;  // total number of steps

    std::vector<Magick::Image> frames;

    int N = 10000;
    std::vector<Body> bodies = generate_random_bodies(N);

    auto start_total = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < steps; ++step) {
        if (algo == 1) {
            sequential_simulation(bodies, dt);
        } else if (algo == 2){
            parallel_step_simulation(bodies, dt);
        }
        else if (algo == 3) {
            parallel_distinc_simulation(bodies, dt);
        }
        else if (algo == 4) {
            parallel_combined_simulation(bodies, dt);
        }
        frames.push_back(drawFrame(bodies));
    }
    
    auto end_total = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = end_total - start_total;
    std::cout << "Total simulation time: " << total_duration.count() << " seconds." << std::endl;

    for (auto& frame : frames) {
        frame.animationDelay(2);
    }
    Magick::writeImages(frames.begin(), frames.end(), "simulation.gif");

    return 0;
}
