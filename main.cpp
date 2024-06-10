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
double THETA = 0.5;

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
            }
        }
    }
}

void compute_forces_not_optimized_naive(std::vector<Body>& bodies, std::vector<double>& fx, std::vector<double>& fy) {
    int n = bodies.size();
    for (int i = 0; i < n; ++i) {
        fx[i] = 0;
        fy[i] = 0;
    }

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

//--------------------------------PARALLEL STEP------------------------------------------------------

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

//--------------------------------PARALLEL FORCES------------------------------------------------------

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

                fx[i] += fx_i;
                fy[i] += fy_i;
            }
        }
    }
}

void parallel_forces_simulation(std::vector<Body>& bodies, double dt, int numThreads) {
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

    for (int i = 0; i < n; ++i) {
        bodies[i].vx += fx[i] / bodies[i].mass * dt;
        bodies[i].vy += fy[i] / bodies[i].mass * dt;
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
}

//--------------------------------PARALLEL COMBINED------------------------------------------------------

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

//--------------------------------GENERATE RANDOM------------------------------------------------------

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

//--------------------------------COLLISION------------------------------------------------------

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

//--------------------------------DRAW FRAME------------------------------------------------------

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

//--------------------------------MAIN------------------------------------------------------

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <1 for sequential_simulation, 2 for parallel_step_simulation, 3 for parallel_distinc_simulation, 4 for parallel_combined_simulation, 5 for barnes_hutt_simulation>, and <0 for no collision and 1 for collision> " << std::endl;
        return 1;
    }

    int algo = std::atoi(argv[1]);
    if (algo != -2 && algo != -1 && algo != 0 && algo != 1 && algo != 2 && algo != 3 && algo != 4 && algo != 5) {
        std::cerr << "Invalid option. Use 1 for sequential_simulation, 2 for parallel_step_simulation, 3 for parallel_distinc_simulation, 4 for parallel_combined_simulation, 5 for barnes_hutt_simulation." << std::endl;
        return 1;
    }

    int collision = std::atoi(argv[2]);
    if (collision != 1 && collision != 0) {
        std::cerr << "Invalid option. Use 1 for collision and 0 for no collision in the simulation." << std::endl;
        return 1;
    }
    if (algo == -2) {

        int N = 50;

        int steps = 1000;
        double dt = 1e7;
        std::vector<Body> original_bodies = generate_random_bodies(N);
        std::vector<double> theta_values = {0.1, 0.3, 0.5, 0.7, 0.9};

        nlohmann::json results = nlohmann::json::array();

        for (double theta : theta_values) {
            std::cout<<theta<<std::endl;
            THETA = theta;
            std::vector<Body> sequential_bodies = original_bodies;
            std::vector<Body> barnes_hutt_bodies = original_bodies;

            nlohmann::json sequential_result;
            sequential_result["theta"] = theta;
            sequential_result["trajectory"] = nlohmann::json::array();

            auto start_sequential = std::chrono::high_resolution_clock::now();
            for (int step = 0; step < steps; ++step) {
                parallel_combined_simulation(sequential_bodies, dt, std::thread::hardware_concurrency());
            }
            for (size_t idx = 0; idx < sequential_bodies.size(); ++idx) {
                sequential_result["trajectory"].push_back({{"x", sequential_bodies[idx].x}, {"y", sequential_bodies[idx].y}});
            }
            auto end_sequential = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> total_duration_sequential = end_sequential - start_sequential;
            sequential_result["algorithm"] = "parallel_combined_simulation";
            sequential_result["time"] = total_duration_sequential.count();
            results.push_back(sequential_result);

            nlohmann::json barnes_hutt_result;
            barnes_hutt_result["theta"] = theta;
            barnes_hutt_result["trajectory"] = nlohmann::json::array();

            auto start_barnes_hutt = std::chrono::high_resolution_clock::now();
            for (int step = 0; step < steps; ++step) {
                QuadNode root(-1e12, -1e12, 2e12, 2e12);
                for (auto& body : barnes_hutt_bodies) {
                    root.insert(&body);
                }
                barnes_hutt_simulation(barnes_hutt_bodies, root, dt);
            }
            for (size_t idx = 0; idx < barnes_hutt_bodies.size(); ++idx) {
                barnes_hutt_result["trajectory"].push_back({{"x", barnes_hutt_bodies[idx].x}, {"y", barnes_hutt_bodies[idx].y}});
            }
            auto end_barnes_hutt = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> total_duration_barnes_hutt = end_barnes_hutt - start_barnes_hutt;
            barnes_hutt_result["algorithm"] = "barnes_hutt_simulation";
            barnes_hutt_result["time"] = total_duration_barnes_hutt.count();
            results.push_back(barnes_hutt_result);

        }

        std::ofstream file("simulation_accuracy_results.json");
        file << results.dump(4);
        file.close();


    }
    else if (algo == 0) {
        int N = 100000;
        double dt = 1e7;  // time step in seconds
        int steps = 1;  // total number of steps
        std::vector<Body> bodies = generate_random_bodies(N);
        nlohmann::json results = nlohmann::json::array();
        for (int numThreads = 1; numThreads <= std::thread::hardware_concurrency(); ++numThreads) {
            std::cout<<"numThreads: " << numThreads <<std::endl;
            for (int algoTest = 1; algoTest <= 5; ++algoTest) {std::vector<Body> bodiesCopy = bodies;
                std::cout<<" algo: " << algoTest <<std::endl;
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
                        parallel_forces_simulation(bodiesCopy, dt, numThreads);
                    } else if (algoTest == 4) {
                        parallel_combined_simulation(bodiesCopy, dt, numThreads);
                    }
                    else if (algoTest == 5) {
                        barnes_hutt_simulation(bodiesCopy, root, dt);
                    }
                    if (collision == 1) {
                        handle_collisions(bodies);
                    }
                }
                auto end_total = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> total_duration = end_total - start_total;

                std::cout<<"numThreads: " << numThreads <<" algo: " << algoTest <<" time" << total_duration.count() <<std::endl;

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
    else if (algo == -1) {
        double dt = 1e7;   // time step in seconds
        int steps = 1000;   // total number of steps

        std::vector<Body> bodies = generate_random_bodies(50);

        nlohmann::json results = nlohmann::json::array();
        for (int algo = 1; algo <= 5; ++algo) {
            std::vector<Body> bodiesCopy = bodies;
            nlohmann::json algoResults;
            algoResults["algorithm"] = algo;

            QuadNode root(-1e12, -1e12, 2e12, 2e12);
            if (algo == 5) {
                for (auto& body : bodiesCopy) {
                    root.insert(&body);
                }
            }

            for (int step = 0; step < steps; ++step) {
                if (algo == 1) {
                    sequential_simulation(bodiesCopy, dt);
                }
                else if (algo == 2) {
                    parallel_step_simulation(bodiesCopy, dt, std::thread::hardware_concurrency());
                }
                else if (algo == 3) {
                    parallel_forces_simulation(bodiesCopy, dt, std::thread::hardware_concurrency());
                }
                else if (algo == 4) {
                    parallel_combined_simulation(bodiesCopy, dt, std::thread::hardware_concurrency());
                }
                else if (algo == 5) {
                    barnes_hutt_simulation(bodiesCopy, root, dt);
                }
                if (collision == 1) {
                    handle_collisions(bodies);
                }
                for (size_t idx = 0; idx < bodiesCopy.size(); ++idx) {
                    algoResults["body" + std::to_string(idx)].push_back({{"x", bodiesCopy[idx].x}, {"y", bodiesCopy[idx].y}});
                }
            }
            results.push_back(algoResults);
        }

        std::ofstream file("trajectory_results.json");
        file << results.dump(4);
        file.close();
    }
    else{
        std::vector<Body> bodies = {
            {1.989e30, 0, 0, 0, 0, "yellow", std::log(1.989e30)/10}, // mass, x, y, vx, vy, random color, radius (we define it as log(mass)/10)
            {3.285e23, 7e10, 0, 0, -42000, "#B7B8B9", std::log(3.285e23)/10},
            {4.867e24, -11e10, 0, 0, 36000, "#928590", std::log(4.867e24)/10},
            {5.972e24, 16e10, 0, 0, -29000, "#6b93d6", std::log(5.972e24)/10},
            {6.39e23, -23e10, 0, 0, 24000, "#c1440e", std::log(5.972e24)/10}
        };

        // std::vector<Body> bodies = {
        //     {1e24, 0, 0, 0, 0, "white", std::log(1e24)/10}, // mass, x, y, vx, vy, color, radius 
        //     {1e2, 5e10, 0, 10, 0, "red", std::log(1e2)/10}
        // };

        // std::vector<Body> bodies = generate_random_bodies(100);

        double dt = 5e4;  // time step in seconds
        int steps = 1000;  // total number of steps

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
                parallel_forces_simulation(bodies, dt, numThreads);
            }
            else if (algo == 4) {
                parallel_combined_simulation(bodies, dt, numThreads);
            }
            else if (algo == 5) {
                barnes_hutt_simulation(bodies, root, dt);
            }
            if (collision == 1) {
                handle_collisions(bodies);
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
}
