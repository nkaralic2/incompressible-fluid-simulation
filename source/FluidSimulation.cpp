
#include "FluidSimulation.h"

FluidSimulation::FluidSimulation(
    int nx,
    int ny,
    double cellSize,
    double dt
)
    : nx(nx),
    ny(ny),
    cellSize(cellSize),
    dt(dt),
    gravity(-9.81),
    grid(nx, ny, cellSize)
{
}

void FluidSimulation::initializeParticles() {

    particles.clear();

    // Pocetna oblast fluida
    for (int j = 2; j < ny / 2; j++) {
        for (int i = 2; i < nx - 2; i++) {

            double x = (i + 0.5) * cellSize;
            double y = (j + 0.5) * cellSize;

            particles.emplace_back(x, y);
        }
    }
}

void FluidSimulation::applyForces() {
    // Primjena gravitacije
}

void FluidSimulation::particleToGrid() {
    // P2G: prenos brzina sa cestica na grid
}

void FluidSimulation::computeDivergence() {
    // Racunanje divergencije
}

void FluidSimulation::solvePressure() {
    // Iterativno rjesavanje jednacine pritiska
    // Jacobijeva metoda
}

void FluidSimulation::projectVelocity() {
    // Korekcija brzina pomocu gradijenta pritiska
}

void FluidSimulation::handleBoundaries() {
    // Obrada granica i prepreka
}

void FluidSimulation::gridToParticle() {
    // G2P: interpolacija brzina sa grida na cestice
}

void FluidSimulation::advectParticles() {
    // Azuriranje pozicija cestica
}

void FluidSimulation::step() {

    applyForces();

    particleToGrid();

    computeDivergence();

    solvePressure();

    projectVelocity();

    handleBoundaries();

    gridToParticle();

    advectParticles();
}

void FluidSimulation::run(int steps) {

    for (int i = 0; i < steps; i++) {
        step();
    }
}

const std::vector<Particle>&
FluidSimulation::getParticles() const {
    return particles;
}
