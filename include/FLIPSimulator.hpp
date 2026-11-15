#ifndef FLIP_SIMULATOR_HPP
#define FLIP_SIMULATOR_HPP

/**
 * @file FLIPSimulator.hpp
 * @brief Main simulator class implementing the PIC/FLIP method for 2D incompressible flow.
 *
 * The simulation loop follows these stages each timestep:
 *   1. particleToGrid()        — Transfer particle velocities onto MAC grid (P2G)
 *   2. updateCellTypes()       — Classify cells as FLUID/AIR/SOLID from particle positions
 *   3. savePreviousVelocities()— Snapshot grid velocities for FLIP delta computation
 *   4. applyGravity(dt)        — Add gravitational acceleration to v-faces
 *   5. applyBoundaryConditions() — Zero velocities at solid boundaries
 *   6. solvePressure(dt)       — Pressure projection (Jacobi) + velocity correction
 *   7. applyBoundaryConditions() — Re-enforce solid boundary velocities
 *   8. gridToParticle()        — Transfer corrected velocities back to particles (G2P)
 *   9. updateParticlePositions(dt) — Advect particles (RK2 midpoint method)
 *  10. handleCollisions()      — Clamp particles to domain, kill wall-penetrating velocities
 */

#include "Particle.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"

#include <vector>
#include <string>

/**
 * @struct SimConfig
 * @brief Holds all configurable simulation parameters.
 */
struct SimConfig {
    // Domain
    int gridNX          = 128;          ///< Grid cells in X
    int gridNY          = 128;          ///< Grid cells in Y
    float domainWidth   = 4.0f;        ///< Physical domain width (meters)
    // domainHeight is derived: domainWidth * gridNY / gridNX

    // Particles
    int particlesPerCell = 4;           ///< Particles per cell (2x2 sub-grid)
    float particleRadius = 0.01f;       ///< Collision clearance from walls

    // Dam-break region (fraction of domain)
    float damLeft       = 0.0f;        ///< Left edge of initial water region (fraction)
    float damRight      = 0.35f;       ///< Right edge (fraction of domainWidth)
    float damBottom     = 0.0f;        ///< Bottom edge (fraction)
    float damTop        = 0.6f;        ///< Top edge (fraction of domainHeight)

    // Physics
    float gravity       = -9.81f;      ///< Gravitational acceleration (m/s²)
    float density       = 1000.0f;     ///< Fluid density (kg/m³)
    float dt            = 1.0f / 60.0f;///< Timestep (seconds)
    float flipAlpha     = 0.95f;       ///< PIC/FLIP blending: 0=pure PIC, 1=pure FLIP

    // Solver
    int pressureIterations = 100;      ///< Jacobi iterations for pressure solve
    float pressureTolerance = 1.0e-4f; ///< Relative Jacobi update tolerance

    // Simulation
    int numSteps        = 500;         ///< Total timesteps to simulate
    int outputInterval  = 10;          ///< Write particle data every N steps (0 = disable)
    std::string outputDir = "output";  ///< Directory for CSV output files
    bool enableOutput   = true;        ///< Master switch for file output

    // Random seed for deterministic initialization
    unsigned int seed   = 42;

    /// Derived: domain height
    float domainHeight() const {
        return domainWidth * static_cast<float>(gridNY) / static_cast<float>(gridNX);
    }

    /// Derived: cell size
    float cellSize() const {
        return domainWidth / static_cast<float>(gridNX);
    }

    /// Validate dimensions and physical parameters before constructing a grid.
    bool validate(std::string& error) const;
};

/**
 * @class FLIPSimulator
 * @brief Orchestrates the PIC/FLIP incompressible fluid simulation.
 */
class FLIPSimulator {
public:
    SimConfig config;
    Grid grid;
    PressureSolver pressureSolver;
    std::vector<Particle> particles;

    int currentStep = 0;

    // Divergence metrics (updated each step for reporting)
    float divBeforeProjection = 0.0f;
    float divAfterProjection  = 0.0f;

    /// Construct simulator with given configuration.
    explicit FLIPSimulator(const SimConfig& cfg);

    /// Initialize grid and particles (STAGE 1).
    void initialize();

    /// Run one simulation timestep.
    void step(float dt);

    /// Run the full simulation for config.numSteps steps.
    void run();

    // ---- Individual stage functions ----

    /// STAGE 2: Transfer particle velocities to grid (P2G).
    void particleToGrid();

    /// Classify cells as FLUID/AIR/SOLID based on particle positions.
    void updateCellTypes();

    /// Fill unweighted air-adjacent faces from nearby particle-supported faces.
    void extrapolateVelocities();

    /// Save current grid velocities for FLIP delta computation.
    void savePreviousVelocities();

    /// STAGE 3: Apply gravitational acceleration to vertical grid velocities.
    void applyGravity(float dt);

    /// Enforce solid boundary conditions on grid velocities.
    void applyBoundaryConditions();

    /// STAGE 4: Pressure projection.
    void solvePressure(float dt);

    /// STAGE 5: Transfer corrected grid velocities back to particles (G2P).
    void gridToParticle();

    /// STAGE 6: Advect particles using RK2 midpoint method.
    void updateParticlePositions(float dt);

    /// STAGE 6: Handle particle-wall collisions.
    void handleCollisions();

    /// Write particle positions to CSV file.
    void writeParticleData(int stepNumber) const;

private:
    /// Initialize the grid structure.
    void initializeGrid();

    /// Initialize particles in the dam-break region.
    void initializeParticles();

    /**
     * @brief Interpolate u-velocity at an arbitrary position using bilinear interpolation.
     * @param px World X coordinate.
     * @param py World Y coordinate.
     * @param uField The u-velocity array to sample from.
     * @return Interpolated u value.
     *
     * u-faces are located at (i*dx, (j+0.5)*dx) for i ∈ [0,NX], j ∈ [0,NY-1].
     */
    float interpolateU(float px, float py, const std::vector<float>& uField) const;

    /**
     * @brief Interpolate v-velocity at an arbitrary position using bilinear interpolation.
     * @param px World X coordinate.
     * @param py World Y coordinate.
     * @param vField The v-velocity array to sample from.
     * @return Interpolated v value.
     *
     * v-faces are located at ((i+0.5)*dx, j*dx) for i ∈ [0,NX-1], j ∈ [0,NY].
     */
    float interpolateV(float px, float py, const std::vector<float>& vField) const;
};

#endif // FLIP_SIMULATOR_HPP
