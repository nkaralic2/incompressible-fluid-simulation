/**
 * @file Benchmark.cpp
 * @brief Implementation of performance measurement and CSV export.
 */

#include "Benchmark.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

void Benchmark::reserve(int steps) {
    stepTimings.reserve(steps);
}

void Benchmark::addStepTiming(const StepTiming& t) {
    stepTimings.push_back(t);
}

void Benchmark::computeTotals() {
    totalP2G_ms       = 0.0;
    totalGravity_ms   = 0.0;
    totalPressure_ms  = 0.0;
    totalG2P_ms       = 0.0;
    totalAdvection_ms = 0.0;
    totalSimTime_ms   = 0.0;

    for (const auto& t : stepTimings) {
        totalP2G_ms       += t.particleToGrid_ms;
        totalGravity_ms   += t.applyGravity_ms;
        totalPressure_ms  += t.pressureSolve_ms;
        totalG2P_ms       += t.gridToParticle_ms;
        totalAdvection_ms += t.advection_ms;
        totalSimTime_ms   += t.total_ms;
    }
}

void Benchmark::printSummary() const {
    std::cout << "\n========== BENCHMARK SUMMARY ==========" << std::endl;
    std::cout << "Grid:               " << gridNX << " x " << gridNY << std::endl;
    std::cout << "Particles:          " << numParticles << std::endl;
    std::cout << "Timesteps:          " << numSteps << std::endl;
    std::cout << "Pressure iters:     " << pressureIterations << std::endl;
    std::cout << "FLIP alpha:         " << flipAlpha << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Total simulation:   " << totalSimTime_ms << " ms" << std::endl;
    std::cout << "  Particle-to-Grid: " << totalP2G_ms << " ms ("
              << (totalSimTime_ms > 0 ? 100.0 * totalP2G_ms / totalSimTime_ms : 0) << "%)" << std::endl;
    std::cout << "  Gravity:          " << totalGravity_ms << " ms ("
              << (totalSimTime_ms > 0 ? 100.0 * totalGravity_ms / totalSimTime_ms : 0) << "%)" << std::endl;
    std::cout << "  Pressure solve:   " << totalPressure_ms << " ms ("
              << (totalSimTime_ms > 0 ? 100.0 * totalPressure_ms / totalSimTime_ms : 0) << "%)" << std::endl;
    std::cout << "  Grid-to-Particle: " << totalG2P_ms << " ms ("
              << (totalSimTime_ms > 0 ? 100.0 * totalG2P_ms / totalSimTime_ms : 0) << "%)" << std::endl;
    std::cout << "  Advection+Coll:   " << totalAdvection_ms << " ms ("
              << (totalSimTime_ms > 0 ? 100.0 * totalAdvection_ms / totalSimTime_ms : 0) << "%)" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    double avgStep = (numSteps > 0) ? totalSimTime_ms / numSteps : 0;
    std::cout << "Avg step time:      " << avgStep << " ms" << std::endl;
    std::cout << "Avg FPS:            " << (avgStep > 0 ? 1000.0 / avgStep : 0) << std::endl;
    std::cout << "========================================" << std::endl;
}

void Benchmark::exportCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[Benchmark] Could not open " << filename << " for writing." << std::endl;
        return;
    }

    file << "step,total_ms,p2g_ms,gravity_ms,pressure_ms,g2p_ms,advection_ms\n";
    file << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < stepTimings.size(); ++i) {
        const auto& t = stepTimings[i];
        file << i << ","
             << t.total_ms << ","
             << t.particleToGrid_ms << ","
             << t.applyGravity_ms << ","
             << t.pressureSolve_ms << ","
             << t.gridToParticle_ms << ","
             << t.advection_ms << "\n";
    }
    std::cout << "[Benchmark] Per-step timings written to " << filename << std::endl;
}

void Benchmark::exportSummaryCSV(const std::string& filename) const {
    // Check if file exists (to decide whether to write header)
    bool writeHeader = true;
    {
        std::ifstream check(filename);
        if (check.good()) {
            writeHeader = false;
        }
    }

    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "[Benchmark] Could not open " << filename << " for writing." << std::endl;
        return;
    }

    if (writeHeader) {
        file << "grid_nx,grid_ny,particles,steps,pressure_iters,flip_alpha,"
             << "total_ms,p2g_ms,gravity_ms,pressure_ms,g2p_ms,advection_ms,avg_step_ms\n";
    }

    double avgStep = (numSteps > 0) ? totalSimTime_ms / numSteps : 0;
    file << std::fixed << std::setprecision(4);
    file << gridNX << "," << gridNY << ","
         << numParticles << "," << numSteps << ","
         << pressureIterations << "," << flipAlpha << ","
         << totalSimTime_ms << ","
         << totalP2G_ms << ","
         << totalGravity_ms << ","
         << totalPressure_ms << ","
         << totalG2P_ms << ","
         << totalAdvection_ms << ","
         << avgStep << "\n";

    std::cout << "[Benchmark] Summary appended to " << filename << std::endl;
}
