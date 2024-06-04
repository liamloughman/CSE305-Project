#include <iostream>
#include <vector>
#include <cmath>
#include <Magick++.h>
#include <thread>
#include <mutex>
#include "json.hpp"
#include <fstream>
#include <sstream>
#include <random>
#include <ctime>

const double G = 6.67430e-11;
const double THETA = 0.5;

struct Body {
    double mass;
    double x, y;
    double vx, vy;
    std::string color;
    double radius;
};

//--------------------------------SEQUENTIAL------------------------------------------------------

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

//--------------------------------PARALLEL------------------------------------------------------

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

void parallel_step_simulation(std::vector<Body>& bodies, double dt, int numThreads) {
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);
    compute_forces(bodies, fx, fy);

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

void compute_forces_parallel(std::vector<Body>& bodies, std::vector<double>& fx, std::vector<double>& fy, int start, int end) {
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

void parallel_distinc_simulation(std::vector<Body>& bodies, double dt, int numThreads) {
    int n = bodies.size();
    std::vector<double> fx(n, 0), fy(n, 0);

    std::vector<std::thread> threads;

    for (int i = 0; i < numThreads; ++i) {
        int start = i * n / numThreads;
        int end = (i + 1) * n / numThreads;
        threads.push_back(std::thread(compute_forces_parallel, std::ref(bodies), std::ref(fx), std::ref(fy), start, end));
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

void parallel_combined_simulation(std::vector<Body>& bodies, double dt, int numThreads) {
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

//--------------------------------BARNES-HUTT------------------------------------------------------

struct QuadNode {
    double mass;
    double centerX, centerY;
    double minX, minY, maxX, maxY;
    Body* body;
    QuadNode* NW;
    QuadNode* NE;
    QuadNode* SW;
    QuadNode* SE;

    QuadNode(double minX, double minY, double maxX, double maxY)
        : mass(0), centerX(0), centerY(0), minX(minX), minY(minY), maxX(maxX), maxY(maxY),
          body(nullptr), NW(nullptr), NE(nullptr), SW(nullptr), SE(nullptr) {}

    ~QuadNode() {
        delete NW;
        delete NE;
        delete SW;
        delete SE;
    }

    bool isLeaf() const {
        return !NW && !NE && !SW && !SE;
    }

    void insert(Body* b);
    void computeMassDistribution();
    void computeForce(Body* b, double& fx, double& fy);
};

void QuadNode::insert(Body* b) {
    if (!body && isLeaf()) {
        body = b;
        centerX = b->x;
        centerY = b->y;
        mass = b->mass;
        return;
    }
    if (isLeaf()) {
        NW = new QuadNode(minX, (minY + maxY) / 2, (minX + maxX) / 2, maxY);
        NE = new QuadNode((minX + maxX) / 2, (minY + maxY) / 2, maxX, maxY);
        SW = new QuadNode(minX, minY, (minX + maxX) / 2, (minY + maxY) / 2);
        SE = new QuadNode((minX + maxX) / 2, minY, maxX, (minY + maxY) / 2);
        if (body) {
            if (body->x <= (minX + maxX) / 2) {
                if (body->y <= (minY + maxY) / 2) {
                    SW->insert(body);
                } else {
                    NW->insert(body);
                }
            } else {
                if (body->y <= (minY + maxY) / 2) {
                    SE->insert(body);
                } else {
                    NE->insert(body);
                }
            }
            body = nullptr;
        }
    }
    if (b->x <= (minX + maxX) / 2) {
        if (b->y <= (minY + maxY) / 2) {
            SW->insert(b);
        } else {
            NW->insert(b);
        }
    } else {
        if (b->y <= (minY + maxY) / 2) {
            SE->insert(b);
        } else {
            NE->insert(b);
        }
    }
    if (b->x == centerX && b->y == centerY) {
        b->x += ((std::rand() % 100) - 50) * 1e-9;
        b->y += ((std::rand() % 100) - 50) * 1e-9;
    }
}

void QuadNode::computeMassDistribution() {
    if (isLeaf()) {
        return;
    }
    NW->computeMassDistribution();
    NE->computeMassDistribution();
    SW->computeMassDistribution();
    SE->computeMassDistribution();
    mass = NW->mass + NE->mass + SW->mass + SE->mass;
    centerX = (NW->centerX * NW->mass + NE->centerX * NE->mass + SW->centerX * SW->mass + SE->centerX * SE->mass) / mass;
    centerY = (NW->centerY * NW->mass + NE->centerY * NE->mass + SW->centerY * SW->mass + SE->centerY * SE->mass) / mass;
}

void QuadNode::computeForce(Body* b, double& fx, double& fy) {
    if (isLeaf() && body == b) {
        return;
    }
    double dx = centerX - b->x;
    double dy = centerY - b->y;
    double distSq = std::max(dx * dx + dy * dy, 1e-6);
    double dist = std::sqrt(distSq);  
    if ((maxX - minX) / dist < THETA || isLeaf()) {
        double force = G * mass * b->mass / distSq;
        fx += force * dx / dist;
        fy += force * dy / dist;
        return;
    }
    if (NW) NW->computeForce(b, fx, fy);
    if (NE) NE->computeForce(b, fx, fy);
    if (SW) SW->computeForce(b, fx, fy);
    if (SE) SE->computeForce(b, fx, fy);
}

void barnes_hutt_simulation(std::vector<Body>& bodies, QuadNode& root, double dt) {
    root.computeMassDistribution();
    for (auto& body : bodies) {
        double fx = 0, fy = 0;
        root.computeForce(&body, fx, fy);
        body.vx += fx / body.mass * dt;
        body.vy += fy / body.mass * dt;
    }
    for (auto& body : bodies) {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
    }
}

std::string generate_random_color() {
    std::ostringstream color;
    color << "#";
    for (int i = 0; i < 6; ++i) {
        int value = std::rand() % 16;
        color << std::hex << value;
    }
    return color.str();
}

std::vector<Body> generate_random_bodies(int N) {
    std::vector<Body> bodies;
    bodies.reserve(N);
    std::mt19937_64 rng(std::time(nullptr));
    std::uniform_real_distribution<double> log_mass_dist(2, 24);
    std::uniform_real_distribution<double> pos_dist(-30e10, 30e10);
    std::uniform_real_distribution<double> vel_dist(-100, 100);
    bodies.push_back({1e25, 0, 0, 0, 0, "white", std::log(1e25)/10}); // central massive body
    for (int i = 1; i < N; ++i) {
        double exponent = log_mass_dist(rng);
        double mass = std::pow(10, exponent);
        double x = pos_dist(rng);
        double y = pos_dist(rng);
        double vx = vel_dist(rng);
        double vy = vel_dist(rng);
        bodies.push_back({mass, x, y, vx, vy, generate_random_color(), std::log(mass)/10});
    }
    return bodies;
}

void handle_collisions(std::vector<Body>& bodies) {
    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double distance = std::sqrt(dx * dx + dy * dy);
            double combinedRadius = bodies[i].radius + bodies[j].radius;

            if (distance < combinedRadius * 1.567e9) {
                double nx = dx / distance;
                double ny = dy / distance;
                double dvx = bodies[j].vx - bodies[i].vx;
                double dvy = bodies[j].vy - bodies[i].vy;
                double relativeVelocity = dvx * nx + dvy * ny;

                if (relativeVelocity < 0) {
                    double impulse = -2 * relativeVelocity / (1 / bodies[i].mass + 1 / bodies[j].mass);
                    bodies[i].vx -= impulse / bodies[i].mass * nx;
                    bodies[i].vy -= impulse / bodies[i].mass * ny;
                    bodies[j].vx += impulse / bodies[j].mass * nx;
                    bodies[j].vy += impulse / bodies[j].mass * ny;
                    double overlap = (combinedRadius * 1.567e9 - distance) / 2.0;
                    bodies[i].x -= overlap * nx;
                    bodies[i].y -= overlap * ny;
                    bodies[j].x += overlap * nx;
                    bodies[j].y += overlap * ny;
                }
            }
        }
    }
}

Magick::Image drawFrame(const std::vector<Body>& bodies) {
    Magick::Image image(Magick::Geometry(800, 600), "black");
    image.type(Magick::TrueColorType);

    for (const Body& body : bodies) {
        double x = body.x / 1e9 + 400;
        double y = body.y / 1e9 + 300;
        image.strokeColor(body.color);
        image.fillColor(body.color);
        image.draw(Magick::DrawableCircle(x, y, x + body.radius, y + body.radius));
    }
    return image;
}


int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <1 for sequential_simulation, 2 for parallel_step_simulation, 3 for parallel_distinc_simulation, 4 for parallel_combined_simulation, >" << std::endl;
        return 1;
    }

    int algo = std::atoi(argv[1]);
    if (algo != 0 && algo != 1 && algo != 2 && algo != 3 && algo != 4 && algo != 5) {
        std::cerr << "Invalid option. Use 1 for sequential_simulation, 2 for parallel_step_simulation, 3 for parallel_distinc_simulation, 4 for parallel_combined_simulation, 5 for barnes_hutt_simulation." << std::endl;
        return 1;
    }

    if (algo == 0) {
        int N = 10000;
        double dt = 1e7;  // time step in seconds
        int steps = 1;  // total number of steps
        std::vector<Body> bodies = generate_random_bodies(N);
        nlohmann::json results = nlohmann::json::array();
        for (int numThreads = 1; numThreads <= std::thread::hardware_concurrency(); ++numThreads) {
            for (int algoTest = 1; algoTest <= 5; ++algoTest) {std::vector<Body> bodiesCopy = bodies;
                auto start_total = std::chrono::high_resolution_clock::now();
                double minX = -1e12, minY = -1e12, maxX = 1e12, maxY = 1e12;
                QuadNode root(minX, minY, maxX, maxY);
                if (algoTest == 5){
                    for (auto& body : bodies) {
                        root.insert(&body);
                    }
                }
                for (int step = 0; step < steps; ++step) {
                    if (algoTest == 1) {
                        sequential_simulation(bodiesCopy, dt);
                    } else if (algoTest == 2) {
                        parallel_step_simulation(bodiesCopy, dt, numThreads);
                    } else if (algoTest == 3) {
                        parallel_distinc_simulation(bodiesCopy, dt, numThreads);
                    } else if (algoTest == 4) {
                        parallel_combined_simulation(bodiesCopy, dt, numThreads);
                    }
                    else if (algoTest == 5) {
                        barnes_hutt_simulation(bodiesCopy, root, dt);
                    }
                }
                auto end_total = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> total_duration = end_total - start_total;

                nlohmann::json result;
                result["algorithm"] = algoTest;
                result["threads"] = numThreads;
                result["time"] = total_duration.count();
                results.push_back(result);
            }
        }
        std::string filename = "results" + std::to_string(N) + ".json";
        std::ofstream file(filename);
        file << results;
        file.close();

        return 0;

    }
    else{
        // std::vector<Body> bodies = {
        //     {1e24, 0, 0, 0, 0}, // mass, x, y, vx, vy
        //     {1e2, 5e10, 1e10, -70, 60},
        //     {1e2, -5e10, 5e10, 10, -20}
        // };

        // std::vector<Body> bodies = {
        //     {1e24, 0, 0, 0, 0, "white", std::log(1e24)/10}, // mass, x, y, vx, vy, color, radius 
        //     {1e2, 5e10, 0, 10, 0, "red", std::log(1e2)/10}
        // };

        std::vector<Body> bodies = generate_random_bodies(10000);

        double dt = 1e7;  // time step in seconds
        int steps = 1;  // total number of steps

        int numThreads = std::thread::hardware_concurrency();

        std::vector<Magick::Image> frames;

        double minX = -1e12, minY = -1e12, maxX = 1e12, maxY = 1e12;
        QuadNode root(minX, minY, maxX, maxY);
        if (algo == 5){
            for (auto& body : bodies) {
                root.insert(&body);
            }
        }

        auto start_total = std::chrono::high_resolution_clock::now();
        for (int step = 0; step < steps; ++step) {
            if (algo == 1) {
                sequential_simulation(bodies, dt);
            } else if (algo == 2){
                parallel_step_simulation(bodies, dt, numThreads);
            }
            else if (algo == 3) {
                parallel_distinc_simulation(bodies, dt, numThreads);
            }
            else if (algo == 4) {
                parallel_combined_simulation(bodies, dt, numThreads);
            }
            else if (algo == 5) {
                barnes_hutt_simulation(bodies, root, dt);
            }
            handle_collisions(bodies);
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
}
