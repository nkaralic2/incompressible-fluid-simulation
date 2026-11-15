
#ifndef FLUID_SIMULATION_H
#define FLUID_SIMULATION_H

#include <vector>
#include "Particle.h"
#include "Grid.h"

class FluidSimulation {
private:
    int nx, ny;
    double cellSize;
    double dt;
    double gravity;

    Grid grid;
    std::vector<Particle> particles;

    void applyForces();

    void particleToGrid();

    void computeDivergence();

    void solvePressure();

    void projectVelocity();

    void handleBoundaries();

    void gridToParticle();

    void advectParticles();

public:
    FluidSimulation(
        int nx,
        int ny,
        double cellSize,
        double dt
    );

    void initializeParticles();

    void step();

    void run(int steps);

    const std::vector<Particle>& getParticles() const;
};

#endif
