/**
 * @file main.cpp
 * @brief Entry point for the PIC/FLIP incompressible fluid simulator.
 *
 * Supports three modes:
 *   1. Default simulation (dam-break, 128x128, 500 steps)
 *   2. Benchmark mode (multiple configurations)
 *   3. Test mode (correctness validation)
 *
 * Usage:
 *   fluid_sim                       — Run default simulation
 *   fluid_sim --benchmark           — Run benchmark suite
 *   fluid_sim --test                — Run correctness tests
 *   fluid_sim --nx 64 --ny 64 --steps 200 --alpha 0.95  — Custom config
 *   fluid_sim --no-output           — Disable file output
 */

#include "FLIPSimulator.hpp"
#include "Benchmark.hpp"

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <array>
#include <filesystem>
#include <fstream>

// ============================================================================
// Benchmarked simulation run
// ============================================================================

Benchmark runBenchmarked(SimConfig& cfg) {
    FLIPSimulator sim(cfg);
    sim.initialize();

    Benchmark bench;
    bench.gridNX = cfg.gridNX;
    bench.gridNY = cfg.gridNY;
    bench.numParticles = static_cast<int>(sim.particles.size());
    bench.numSteps = cfg.numSteps;
    bench.pressureIterations = cfg.pressureIterations;
    bench.flipAlpha = cfg.flipAlpha;
    bench.reserve(cfg.numSteps);

    for (int s = 0; s < cfg.numSteps; ++s) {
        sim.currentStep = s;
        StepTiming timing;

        auto t0 = Benchmark::now();

        // P2G
        auto tP2G_start = Benchmark::now();
        sim.particleToGrid();
        sim.updateCellTypes();
        sim.savePreviousVelocities();
        auto tP2G_end = Benchmark::now();
        timing.particleToGrid_ms = Benchmark::elapsedMs(tP2G_start, tP2G_end);

        // Gravity
        auto tGrav_start = Benchmark::now();
        sim.applyGravity(cfg.dt);
        sim.applyBoundaryConditions();
        auto tGrav_end = Benchmark::now();
        timing.applyGravity_ms = Benchmark::elapsedMs(tGrav_start, tGrav_end);

        // Pressure
        auto tPress_start = Benchmark::now();
        sim.solvePressure(cfg.dt);
        sim.applyBoundaryConditions();
        auto tPress_end = Benchmark::now();
        timing.pressureSolve_ms = Benchmark::elapsedMs(tPress_start, tPress_end);

        // G2P
        auto tG2P_start = Benchmark::now();
        sim.gridToParticle();
        auto tG2P_end = Benchmark::now();
        timing.gridToParticle_ms = Benchmark::elapsedMs(tG2P_start, tG2P_end);

        // Advection + Collision
        auto tAdv_start = Benchmark::now();
        sim.updateParticlePositions(cfg.dt);
        sim.handleCollisions();
        auto tAdv_end = Benchmark::now();
        timing.advection_ms = Benchmark::elapsedMs(tAdv_start, tAdv_end);

        auto t1 = Benchmark::now();
        timing.total_ms = Benchmark::elapsedMs(t0, t1);

        bench.addStepTiming(timing);

        // Progress (don't print too often for benchmarks)
        if ((s + 1) % 100 == 0 || s == 0 || s == cfg.numSteps - 1) {
            std::cout << "[Step " << (s + 1) << "/" << cfg.numSteps
                      << "] divBefore=" << std::scientific << sim.divBeforeProjection
                      << "  divAfter=" << sim.divAfterProjection
                      << std::fixed << std::endl;
        }
    }

    bench.computeTotals();
    return bench;
}

// ============================================================================
// Correctness tests
// ============================================================================

bool testParticlesBoundary(const SimConfig& cfg) {
    std::cout << "\n--- Test: Particles stay inside domain ---" << std::endl;
    SimConfig testCfg = cfg;
    testCfg.gridNX = 64;
    testCfg.gridNY = 64;
    testCfg.numSteps = 100;
    testCfg.enableOutput = false;

    FLIPSimulator sim(testCfg);
    sim.initialize();
    int initialCount = static_cast<int>(sim.particles.size());

    for (int s = 0; s < testCfg.numSteps; ++s) {
        sim.step(testCfg.dt);
    }

    float dx = testCfg.cellSize();
    float domW = testCfg.domainWidth;
    float domH = testCfg.domainHeight();
    bool allInside = true;
    bool noNaN = true;

    for (const auto& p : sim.particles) {
        float margin = dx + testCfg.particleRadius;
        if (p.x < margin || p.x > domW - margin ||
            p.y < margin || p.y > domH - margin) {
            allInside = false;
        }
        if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.vx) || std::isnan(p.vy)) {
            noNaN = false;
        }
        if (std::isinf(p.x) || std::isinf(p.y) || std::isinf(p.vx) || std::isinf(p.vy)) {
            noNaN = false;
        }
    }

    int finalCount = static_cast<int>(sim.particles.size());
    bool countPreserved = (finalCount == initialCount);

    std::cout << "  Particle count preserved: " << (countPreserved ? "PASS" : "FAIL")
              << " (" << initialCount << " -> " << finalCount << ")" << std::endl;
    std::cout << "  All particles inside solid margin: " << (allInside ? "PASS" : "FAIL") << std::endl;
    std::cout << "  No NaN/Inf values: " << (noNaN ? "PASS" : "FAIL") << std::endl;

    return countPreserved && allInside && noNaN;
}

bool testParticleDistribution(const SimConfig& cfg) {
    std::cout << "\n--- Test: Particle distribution and volume proxy ---" << std::endl;
    SimConfig testCfg = cfg;
    testCfg.gridNX = 64;
    testCfg.gridNY = 64;
    testCfg.numSteps = 100;
    testCfg.enableOutput = false;

    FLIPSimulator sim(testCfg);
    sim.initialize();
    const int initialCount = static_cast<int>(sim.particles.size());

    auto measure = [&]() {
        std::vector<int> counts(sim.grid.cellSize(), 0);
        for (const auto& p : sim.particles) {
            int i = std::max(0, std::min(static_cast<int>(p.x / sim.grid.dx), sim.grid.NX - 1));
            int j = std::max(0, std::min(static_cast<int>(p.y / sim.grid.dx), sim.grid.NY - 1));
            ++counts[sim.grid.cellIndex(i, j)];
        }

        int occupied = 0;
        int maximum = 0;
        double sum = 0.0;
        for (int count : counts) {
            if (count > 0) {
                ++occupied;
                maximum = std::max(maximum, count);
                sum += count;
            }
        }
        double mean = occupied > 0 ? sum / occupied : 0.0;
        double variance = 0.0;
        for (int count : counts) {
            if (count > 0) {
                double delta = count - mean;
                variance += delta * delta;
            }
        }
        double coefficientOfVariation = occupied > 0 && mean > 0.0
            ? std::sqrt(variance / occupied) / mean : 0.0;
        return std::array<double, 4>{static_cast<double>(occupied), mean,
                                     static_cast<double>(maximum), coefficientOfVariation};
    };

    const auto initial = measure();
    for (int s = 0; s < testCfg.numSteps; ++s) {
        sim.step(testCfg.dt);
    }
    const auto final = measure();

    std::cout << "  Occupied cells: " << initial[0] << " -> " << final[0] << std::endl;
    std::cout << "  Mean particles/occupied cell: " << initial[1] << " -> " << final[1] << std::endl;
    std::cout << "  Maximum particles/cell: " << initial[2] << " -> " << final[2] << std::endl;
    std::cout << "  Occupancy CV: " << initial[3] << " -> " << final[3] << std::endl;

    bool finite = true;
    for (double value : final) finite = finite && std::isfinite(value);
    bool passed = finite && initialCount == static_cast<int>(sim.particles.size()) && final[0] > 0.0;
    std::cout << "  Particle count and occupancy remain valid: "
              << (passed ? "PASS" : "FAIL") << std::endl;
    return passed;
}

bool testResolutionStability(const SimConfig& cfg) {
    std::cout << "\n--- Test: Multi-resolution numerical stability ---" << std::endl;
    bool passed = true;
    for (int resolution : {32, 64, 128}) {
        SimConfig testCfg = cfg;
        testCfg.gridNX = resolution;
        testCfg.gridNY = resolution;
        testCfg.numSteps = 30;
        testCfg.enableOutput = false;

        FLIPSimulator sim(testCfg);
        sim.initialize();
        const int initialCount = static_cast<int>(sim.particles.size());
        float peakAfter = 0.0f;
        for (int s = 0; s < testCfg.numSteps; ++s) {
            sim.step(testCfg.dt);
            peakAfter = std::max(peakAfter, sim.divAfterProjection);
        }

        bool finite = std::isfinite(peakAfter);
        bool countPreserved = initialCount == static_cast<int>(sim.particles.size());
        bool inside = true;
        float margin = testCfg.cellSize() + testCfg.particleRadius;
        for (const auto& p : sim.particles) {
            inside = inside && p.x >= margin && p.x <= testCfg.domainWidth - margin &&
                     p.y >= margin && p.y <= testCfg.domainHeight() - margin;
        }
        bool resolutionPassed = finite && countPreserved && inside;
        passed = passed && resolutionPassed;
        std::cout << "  " << resolution << "x" << resolution
                  << ": peak divAfter=" << peakAfter
                  << ", count=" << (countPreserved ? "PASS" : "FAIL")
                  << ", boundary=" << (inside ? "PASS" : "FAIL") << std::endl;
    }
    return passed;
}

bool testOutputFrameLifecycle(const SimConfig& cfg) {
    std::cout << "\n--- Test: Output frame lifecycle ---" << std::endl;
    namespace fs = std::filesystem;
    SimConfig testCfg = cfg;
    testCfg.gridNX = 16;
    testCfg.gridNY = 16;
    testCfg.numSteps = 8;
    testCfg.outputInterval = 2;
    testCfg.outputDir = "test_output_frames";
    testCfg.enableOutput = true;

    fs::create_directories(testCfg.outputDir);
    std::ofstream stale(testCfg.outputDir + "/particles_99999.csv");
    stale << "x,y,vx,vy\n";
    stale.close();

    FLIPSimulator sim(testCfg);
    sim.initialize();
    sim.run();

    int frameCount = 0;
    for (const fs::directory_entry& entry : fs::directory_iterator(testCfg.outputDir)) {
        const std::string name = entry.path().filename().string();
        if (entry.is_regular_file() && name.rfind("particles_", 0) == 0 &&
            entry.path().extension() == ".csv") {
            ++frameCount;
        }
    }
    fs::remove_all(testCfg.outputDir);

    const int expectedFrames = 5; // steps 0, 2, 4, 6, and the final step 7
    bool passed = frameCount == expectedFrames;
    std::cout << "  Frames: " << frameCount << " (expected " << expectedFrames << "): "
              << (passed ? "PASS" : "FAIL") << std::endl;
    return passed;
}

bool testGravity(const SimConfig& cfg) {
    std::cout << "\n--- Test: Gravity moves water downward ---" << std::endl;
    SimConfig testCfg = cfg;
    testCfg.gridNX = 32;
    testCfg.gridNY = 32;
    testCfg.numSteps = 50;
    testCfg.enableOutput = false;
    testCfg.damBottom = 0.3f;
    testCfg.damTop = 0.7f;
    testCfg.damLeft = 0.3f;
    testCfg.damRight = 0.7f;

    FLIPSimulator sim(testCfg);
    sim.initialize();

    // Compute initial center of mass
    float initialComy = 0.0f;
    for (const auto& p : sim.particles) {
        initialComy += p.y;
    }
    initialComy /= static_cast<float>(sim.particles.size());

    for (int s = 0; s < testCfg.numSteps; ++s) {
        sim.step(testCfg.dt);
    }

    float finalComy = 0.0f;
    for (const auto& p : sim.particles) {
        finalComy += p.y;
    }
    finalComy /= static_cast<float>(sim.particles.size());

    bool moved_down = (finalComy < initialComy);
    std::cout << "  Initial CoM_y: " << initialComy << std::endl;
    std::cout << "  Final CoM_y:   " << finalComy << std::endl;
    std::cout << "  Water moved downward: " << (moved_down ? "PASS" : "FAIL") << std::endl;

    return moved_down;
}

bool testDivergenceReduction(const SimConfig& cfg) {
    std::cout << "\n--- Test: Pressure reduces divergence ---" << std::endl;
    SimConfig testCfg = cfg;
    testCfg.gridNX = 64;
    testCfg.gridNY = 64;
    testCfg.numSteps = 10;
    testCfg.enableOutput = false;

    FLIPSimulator sim(testCfg);
    sim.initialize();

    bool allReduced = true;
    for (int s = 0; s < testCfg.numSteps; ++s) {
        sim.step(testCfg.dt);
        // After step 0 (which may start from zero divergence), check later steps
        if (s > 0 && sim.divBeforeProjection > 1e-6f) {
            if (sim.divAfterProjection >= sim.divBeforeProjection) {
                std::cout << "  Step " << s << ": divergence NOT reduced ("
                          << sim.divBeforeProjection << " -> " << sim.divAfterProjection << ")" << std::endl;
                allReduced = false;
            }
        }
    }

    std::cout << "  Last step: divBefore=" << sim.divBeforeProjection
              << " divAfter=" << sim.divAfterProjection << std::endl;
    std::cout << "  Pressure iterations: " << sim.pressureSolver.lastIterations
              << ", pressure delta=" << sim.pressureSolver.lastPressureDelta
              << ", Poisson residual=" << sim.pressureSolver.lastResidual
              << " (ratio=" << sim.pressureSolver.residualRatio << ")" << std::endl;
    std::cout << "  Pressure reduces divergence: " << (allReduced ? "PASS" : "FAIL") << std::endl;

    return allReduced;
}

bool testDeterminism(const SimConfig& cfg) {
    std::cout << "\n--- Test: Deterministic results ---" << std::endl;
    SimConfig testCfg = cfg;
    testCfg.gridNX = 32;
    testCfg.gridNY = 32;
    testCfg.numSteps = 20;
    testCfg.enableOutput = false;

    // Run 1
    FLIPSimulator sim1(testCfg);
    sim1.initialize();
    for (int s = 0; s < testCfg.numSteps; ++s) sim1.step(testCfg.dt);

    // Run 2
    FLIPSimulator sim2(testCfg);
    sim2.initialize();
    for (int s = 0; s < testCfg.numSteps; ++s) sim2.step(testCfg.dt);

    bool identical = (sim1.particles.size() == sim2.particles.size());
    if (identical) {
        for (size_t i = 0; i < sim1.particles.size(); ++i) {
            if (sim1.particles[i].x != sim2.particles[i].x ||
                sim1.particles[i].y != sim2.particles[i].y ||
                sim1.particles[i].vx != sim2.particles[i].vx ||
                sim1.particles[i].vy != sim2.particles[i].vy) {
                identical = false;
                break;
            }
        }
    }

    std::cout << "  Results identical across runs: " << (identical ? "PASS" : "FAIL") << std::endl;
    return identical;
}

bool testSolidBoundaryVelocities(const SimConfig& cfg) {
    std::cout << "\n--- Test: Solid boundary normal velocities ---" << std::endl;
    SimConfig testCfg = cfg;
    testCfg.gridNX = 32;
    testCfg.gridNY = 32;
    testCfg.numSteps = 20;
    testCfg.enableOutput = false;

    FLIPSimulator sim(testCfg);
    sim.initialize();
    for (int s = 0; s < testCfg.numSteps; ++s) {
        sim.step(testCfg.dt);
    }

    bool boundaryZero = true;
    for (int j = 0; j < sim.grid.NY; ++j) {
        boundaryZero = boundaryZero && sim.grid.u[sim.grid.uIndex(0, j)] == 0.0f;
        boundaryZero = boundaryZero && sim.grid.u[sim.grid.uIndex(sim.grid.NX, j)] == 0.0f;
    }
    for (int i = 0; i < sim.grid.NX; ++i) {
        boundaryZero = boundaryZero && sim.grid.v[sim.grid.vIndex(i, 0)] == 0.0f;
        boundaryZero = boundaryZero && sim.grid.v[sim.grid.vIndex(i, sim.grid.NY)] == 0.0f;
    }

    std::cout << "  Solid boundary normals are zero: "
              << (boundaryZero ? "PASS" : "FAIL") << std::endl;
    return boundaryZero;
}

bool testPICFLIPBlending(const SimConfig& cfg) {
    std::cout << "\n--- Test: PIC/FLIP blending values ---" << std::endl;
    bool allOk = true;

    for (float alpha : {0.0f, 0.5f, 0.95f, 1.0f}) {
        SimConfig testCfg = cfg;
        testCfg.gridNX = 32;
        testCfg.gridNY = 32;
        testCfg.numSteps = 10;
        testCfg.enableOutput = false;
        testCfg.flipAlpha = alpha;

        FLIPSimulator sim(testCfg);
        sim.initialize();

        bool crashed = false;
        for (int s = 0; s < testCfg.numSteps; ++s) {
            sim.step(testCfg.dt);
            for (const auto& p : sim.particles) {
                if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.vx) || std::isnan(p.vy)) {
                    crashed = true;
                    break;
                }
            }
            if (crashed) break;
        }

        std::cout << "  alpha=" << alpha << ": " << (crashed ? "FAIL (NaN)" : "PASS") << std::endl;
        if (crashed) allOk = false;
    }

    return allOk;
}

bool testPressureSolverSmallGrid() {
    std::cout << "\n--- Test: Pressure solver on small grid ---" << std::endl;

    // Create a small 4x4 grid with known conditions
    // Set up a single fluid cell with known divergence and check that
    // the pressure solver produces a reasonable result
    Grid grid(6, 6, 1.0f);

    // Set boundary as solid
    for (int j = 0; j < 6; ++j) {
        for (int i = 0; i < 6; ++i) {
            if (i == 0 || i == 5 || j == 0 || j == 5) {
                grid.cellType[grid.cellIndex(i, j)] = CellType::SOLID;
            } else {
                grid.cellType[grid.cellIndex(i, j)] = CellType::FLUID;
            }
        }
    }

    // Set up a velocity field with known divergence
    // u[i+1,j] - u[i,j] + v[i,j+1] - v[i,j] = divergence * dx
    // Create a divergent field in cell (3,3): u_right > u_left
    grid.u[grid.uIndex(4, 3)] = 1.0f;  // right face of cell (3,3)
    grid.u[grid.uIndex(3, 3)] = 0.0f;  // left face
    grid.v[grid.vIndex(3, 4)] = 0.0f;  // top face
    grid.v[grid.vIndex(3, 3)] = 0.0f;  // bottom face

    PressureSolver solver(200, 1000.0f);

    // Compute initial divergence
    solver.computeDivergence(grid);
    float divBefore = solver.computeMaxDivergence(grid);

    // Solve
    float dt = 1.0f / 60.0f;
    solver.solve(grid, dt);

    // Recompute divergence
    solver.computeDivergence(grid);
    float divAfter = solver.computeMaxDivergence(grid);

    std::cout << "  Divergence before: " << divBefore << std::endl;
    std::cout << "  Divergence after:  " << divAfter << std::endl;
    std::cout << "  Pressure iterations: " << solver.lastIterations
              << ", pressure delta=" << solver.lastPressureDelta
              << ", Poisson residual=" << solver.lastResidual
              << " (ratio=" << solver.residualRatio << ")" << std::endl;

    bool reduced = (divAfter < divBefore) || (divBefore < 1e-10f);
    std::cout << "  Divergence reduced: " << (reduced ? "PASS" : "FAIL") << std::endl;

    // Check that pressure is not all zero
    float maxP = 0.0f;
    for (int j = 1; j < 5; ++j) {
        for (int i = 1; i < 5; ++i) {
            float p = std::fabs(grid.pressure[grid.cellIndex(i, j)]);
            if (p > maxP) maxP = p;
        }
    }
    bool hasNonzero = maxP > 1e-10f;
    std::cout << "  Max pressure: " << maxP << " (non-zero: " << (hasNonzero ? "PASS" : "FAIL") << ")" << std::endl;

    return reduced && hasNonzero;
}

bool testConfigurationValidation() {
    std::cout << "\n--- Test: Configuration validation ---" << std::endl;
    SimConfig invalid;
    invalid.gridNX = 2;

    std::string error;
    bool rejectedByValidator = !invalid.validate(error) && !error.empty();
    bool rejectedBySimulator = false;
    try {
        FLIPSimulator sim(invalid);
    } catch (const std::invalid_argument&) {
        rejectedBySimulator = true;
    }

    bool passed = rejectedByValidator && rejectedBySimulator;
    std::cout << "  Invalid configuration rejected: "
              << (passed ? "PASS" : "FAIL") << std::endl;
    return passed;
}

int runTests() {
    std::cout << "========================================" << std::endl;
    std::cout << "   FLIP Simulator Correctness Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    SimConfig baseCfg;
    int passed = 0, total = 0;

    auto check = [&](bool result) { total++; if (result) passed++; };

    check(testParticlesBoundary(baseCfg));
    check(testParticleDistribution(baseCfg));
    check(testGravity(baseCfg));
    check(testDivergenceReduction(baseCfg));
    check(testSolidBoundaryVelocities(baseCfg));
    check(testResolutionStability(baseCfg));
    check(testDeterminism(baseCfg));
    check(testPICFLIPBlending(baseCfg));
    check(testPressureSolverSmallGrid());
    check(testConfigurationValidation());
    check(testOutputFrameLifecycle(baseCfg));

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Results: " << passed << "/" << total << " tests passed" << std::endl;
    std::cout << "========================================" << std::endl;

    return (passed == total) ? 0 : 1;
}

// ============================================================================
// Benchmark suite
// ============================================================================

int runBenchmarks() {
    std::cout << "========================================" << std::endl;
    std::cout << "   FLIP Simulator Benchmark Suite" << std::endl;
    std::cout << "========================================" << std::endl;

    struct BenchConfig {
        int nx, ny, steps;
        std::string label;
    };

    std::vector<BenchConfig> configs = {
        {64,  64,  500, "64x64"},
        {128, 128, 500, "128x128"},
        {256, 256, 500, "256x256"},
    };

    for (const auto& bc : configs) {
        std::cout << "\n--- Benchmark: " << bc.label << " ---" << std::endl;

        SimConfig cfg;
        cfg.gridNX = bc.nx;
        cfg.gridNY = bc.ny;
        cfg.numSteps = bc.steps;
        cfg.enableOutput = false;

        Benchmark bench = runBenchmarked(cfg);
        bench.printSummary();
        bench.exportCSV("benchmarks/bench_" + bc.label + ".csv");
        bench.exportSummaryCSV("benchmarks/benchmark_summary.csv");
    }

    return 0;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    SimConfig cfg;
    bool benchmarkMode = false;
    bool testMode = false;
    bool doBenchmarkTiming = false;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--benchmark") {
            benchmarkMode = true;
        } else if (arg == "--test") {
            testMode = true;
        } else if (arg == "--timed") {
            doBenchmarkTiming = true;
        } else if (arg == "--no-output") {
            cfg.enableOutput = false;
        } else if (arg == "--nx" && i + 1 < argc) {
            cfg.gridNX = std::stoi(argv[++i]);
        } else if (arg == "--ny" && i + 1 < argc) {
            cfg.gridNY = std::stoi(argv[++i]);
        } else if (arg == "--steps" && i + 1 < argc) {
            cfg.numSteps = std::stoi(argv[++i]);
        } else if (arg == "--alpha" && i + 1 < argc) {
            cfg.flipAlpha = std::stof(argv[++i]);
        } else if (arg == "--dt" && i + 1 < argc) {
            cfg.dt = std::stof(argv[++i]);
        } else if (arg == "--pressure-iters" && i + 1 < argc) {
            cfg.pressureIterations = std::stoi(argv[++i]);
        } else if (arg == "--pressure-tolerance" && i + 1 < argc) {
            cfg.pressureTolerance = std::stof(argv[++i]);
        } else if (arg == "--output-interval" && i + 1 < argc) {
            cfg.outputInterval = std::stoi(argv[++i]);
        } else if (arg == "--output-dir" && i + 1 < argc) {
            cfg.outputDir = argv[++i];
        } else if (arg == "--seed" && i + 1 < argc) {
            cfg.seed = static_cast<unsigned int>(std::stoi(argv[++i]));
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: fluid_sim [OPTIONS]\n"
                      << "\nOptions:\n"
                      << "  --benchmark         Run benchmark suite\n"
                      << "  --test              Run correctness tests\n"
                      << "  --timed             Run with per-step timing\n"
                      << "  --no-output         Disable file output\n"
                      << "  --nx N              Grid cells in X (default: 128)\n"
                      << "  --ny N              Grid cells in Y (default: 128)\n"
                      << "  --steps N           Number of timesteps (default: 500)\n"
                      << "  --alpha F           PIC/FLIP blend (0=PIC, 1=FLIP, default: 0.95)\n"
                      << "  --dt F              Timestep in seconds (default: 1/60)\n"
                      << "  --pressure-iters N  Jacobi iterations (default: 100)\n"
                      << "  --pressure-tolerance F  Relative pressure convergence tolerance\n"
                      << "  --output-interval N Write every N steps (default: 10)\n"
                      << "  --output-dir DIR    Output directory (default: output)\n"
                      << "  --seed N            Random seed (default: 42)\n"
                      << "  -h, --help          Show this help\n";
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            return 1;
        }
    }

    std::string configError;
    if (!cfg.validate(configError)) {
        std::cerr << "Invalid configuration: " << configError << std::endl;
        return 2;
    }

    if (testMode) {
        return runTests();
    }

    if (benchmarkMode) {
        return runBenchmarks();
    }

    // Default or timed simulation
    if (doBenchmarkTiming) {
        Benchmark bench = runBenchmarked(cfg);
        bench.printSummary();
        bench.exportCSV("benchmarks/bench_custom.csv");
        bench.exportSummaryCSV("benchmarks/benchmark_summary.csv");
    } else {
        FLIPSimulator sim(cfg);
        sim.initialize();
        sim.run();

        std::cout << "\nSimulation complete. " << sim.particles.size() << " particles, "
                  << cfg.numSteps << " steps." << std::endl;
    }

    return 0;
}
