#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

/**
 * @file Benchmark.hpp
 * @brief Performance measurement utilities for the FLIP simulator.
 *
 * Uses std::chrono::steady_clock for timing individual simulation stages
 * and records results in CSV format.
 */

#include <chrono>
#include <string>
#include <vector>

/**
 * @struct StepTiming
 * @brief Timing data for a single simulation step, broken down by stage.
 */
struct StepTiming {
    double particleToGrid_ms = 0.0;
    double applyGravity_ms   = 0.0;
    double pressureSolve_ms  = 0.0;
    double gridToParticle_ms = 0.0;
    double advection_ms      = 0.0;
    double total_ms          = 0.0;
};

/**
 * @class Benchmark
 * @brief Collects and exports per-step timing and aggregate statistics.
 */
class Benchmark {
public:
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;

    std::vector<StepTiming> stepTimings;

    // Aggregate totals (ms)
    double totalSimTime_ms       = 0.0;
    double totalP2G_ms           = 0.0;
    double totalGravity_ms       = 0.0;
    double totalPressure_ms      = 0.0;
    double totalG2P_ms           = 0.0;
    double totalAdvection_ms     = 0.0;

    // Configuration metadata
    int gridNX = 0, gridNY = 0;
    int numParticles = 0;
    int numSteps = 0;
    int pressureIterations = 0;
    float flipAlpha = 0.0f;

    /// Reserve space for step timings.
    void reserve(int steps);

    /// Record timing for one step.
    void addStepTiming(const StepTiming& t);

    /// Compute aggregate totals from step timings.
    void computeTotals();

    /// Print a summary to stdout.
    void printSummary() const;

    /// Export detailed per-step timings to CSV.
    void exportCSV(const std::string& filename) const;

    /// Export a summary row (for comparing configurations) to CSV.
    void exportSummaryCSV(const std::string& filename) const;

    /// Helper: get current time point.
    static TimePoint now() { return Clock::now(); }

    /// Helper: compute elapsed milliseconds between two time points.
    static double elapsedMs(TimePoint start, TimePoint end) {
        return std::chrono::duration<double, std::milli>(end - start).count();
    }
};

#endif // BENCHMARK_HPP
