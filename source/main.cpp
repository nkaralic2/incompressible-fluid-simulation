
#include <iostream>
#include <chrono>

#include "FluidSimulation.h"

int main() {

    // Dimenzije grida
    const int nx = 128;
    const int ny = 128;

    // Velicina jedne grid celije
    const double cellSize = 1.0 / 128.0;

    // Vremenski korak
    const double dt = 1.0 / 60.0;

    // Broj simulacijskih koraka
    const int steps = 100;

    // Kreiranje simulacije
    FluidSimulation simulation(
        nx,
        ny,
        cellSize,
        dt
    );

    // Inicijalizacija cestica
    simulation.initializeParticles();

    std::cout << "Fluid simulation started.\n";

    // Mjerenje vremena
    auto start = std::chrono::steady_clock::now();

    simulation.run(steps);

    auto end = std::chrono::steady_clock::now();

    // Izracunavanje vremena u milisekundama
    double elapsed =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();

    std::cout << "Simulation finished.\n";

    std::cout << "Particles: "
        << simulation.getParticles().size()
        << '\n';

    std::cout << "Simulation steps: "
        << steps << '\n';

    std::cout << "Execution time: "
        << elapsed << " ms\n";

    return 0;
}
