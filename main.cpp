#include <iostream>
#include <vector>
#include <cmath>
#include <Magick++.h>

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
        {20e10, 0, 0, 0, 0}, // mass, x, y, vx, vy
        {1e10, 5e10, 0, 0, 0}
    };

    double dt = 1e7;  // time step in seconds
    int steps = 200;  // total number of steps

    std::vector<std::string> frameFiles;
    for (int step = 0; step < steps; ++step) {
        updateForcesAndPositions(bodies, dt);
        drawFrame(bodies, step);
        frameFiles.push_back("frame" + std::to_string(step) + ".png");
    }

    std::vector<Magick::Image> frames;
    for (const auto& frameFile : frameFiles) {
        Magick::Image frame(frameFile);
        frame.animationDelay(2);
        frames.push_back(Magick::Image(frameFile));
    }
    Magick::writeImages(frames.begin(), frames.end(), "simulation.gif");
    
    for (const auto& frameFile : frameFiles) {
        remove(frameFile.c_str());
    }
    return 0;
}