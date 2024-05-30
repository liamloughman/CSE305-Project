#include <iostream>
#include <vector>
#include <cmath>
#include <Magick++.h>
#include <thread>

const double G = 6.67430e-11;

struct Body {
    double mass;
    double x, y;
    double vx, vy;
};

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

void updateVelocities(std::vector<Body>& bodies, double dt, std::vector<double> fx, std::vector<double> fy, int start, int end) {
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

void drawFrame(const std::vector<Body>& bodies, int frameNumber) {
    Magick::Image image(Magick::Geometry(800, 600), "black");
    image.type(Magick::TrueColorType);
    image.strokeColor("white");
    image.fillColor("white");
    for (const Body& body : bodies) {
        double x = body.x / 1e9 + 400;
        double y = body.y / 1e9 + 300;
        image.draw(Magick::DrawableCircle(x, y, x + 2, y + 2));
    }
    image.write("frame" + std::to_string(frameNumber) + ".png");
}

int main() {
    
    std::vector<Body> bodies = {
        {1e25, 0, 0, 0, 0}, // mass, x, y, vx, vy
        {5e2, 5e10, 0, -150, 80}
    };

    double dt = 5e7;  // time step in seconds
    int steps = 200;  // total number of steps

    std::vector<std::string> frameFiles;
    for (int step = 0; step < steps; ++step) {
        // Seq(bodies, dt);
        Step_Parallel(bodies, dt);
        drawFrame(bodies, step);
        frameFiles.push_back("frame" + std::to_string(step) + ".png");
    }

    std::vector<Magick::Image> frames;
    for (const auto& frameFile : frameFiles) {
        Magick::Image frame(frameFile);
        frames.push_back(Magick::Image(frameFile));
    }
    Magick::writeImages(frames.begin(), frames.end(), "simulation.gif");
    
    for (const auto& frameFile : frameFiles) {
        remove(frameFile.c_str());
    }
    return 0;
}