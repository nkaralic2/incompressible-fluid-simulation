/**
 * @file FLIPSimulator.cpp
 * @brief Implementation of the PIC/FLIP fluid simulator.
 *
 * This file contains the sequential baseline implementation of all six
 * simulation stages: initialization, P2G, external forces, pressure
 * projection, G2P, and particle advection with collision handling.
 */

#include "FLIPSimulator.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <random>
#include <stdexcept>
#include <filesystem>

namespace {

namespace fs = std::filesystem;

SimConfig checkedConfig(const SimConfig& config) {
    std::string error;
    if (!config.validate(error)) {
        throw std::invalid_argument("Invalid simulation configuration: " + error);
    }
    return config;
}

void prepareParticleOutput(const std::string& directory) {
    fs::create_directories(directory);
    for (const fs::directory_entry& entry : fs::directory_iterator(directory)) {
        const std::string filename = entry.path().filename().string();
        if (entry.is_regular_file() &&
            filename.rfind("particles_", 0) == 0 &&
            entry.path().extension() == ".csv") {
            fs::remove(entry.path());
        }
    }
}

}

bool SimConfig::validate(std::string& error) const {
    if (gridNX < 4 || gridNY < 4) {
        error = "gridNX and gridNY must both be at least 4";
    } else if (!(domainWidth > 0.0f) || !std::isfinite(domainWidth)) {
        error = "domainWidth must be finite and positive";
    } else if (particlesPerCell < 1) {
        error = "particlesPerCell must be positive";
    } else if (!(particleRadius >= 0.0f && particleRadius < cellSize())) {
        error = "particleRadius must be non-negative and smaller than cell size";
    } else if (!(dt > 0.0f) || !std::isfinite(dt)) {
        error = "dt must be finite and positive";
    } else if (!(density > 0.0f) || !std::isfinite(density)) {
        error = "density must be finite and positive";
    } else if (pressureIterations < 1) {
        error = "pressureIterations must be positive";
    } else if (!(pressureTolerance > 0.0f) || !std::isfinite(pressureTolerance)) {
        error = "pressureTolerance must be finite and positive";
    } else if (!(flipAlpha >= 0.0f && flipAlpha <= 1.0f)) {
        error = "flipAlpha must be between 0 and 1";
    } else if (numSteps < 0 || outputInterval < 0) {
        error = "numSteps and outputInterval must be non-negative";
    } else if (!(damLeft >= 0.0f && damLeft < damRight && damRight <= 1.0f &&
                 damBottom >= 0.0f && damBottom < damTop && damTop <= 1.0f)) {
        error = "dam bounds must be ordered fractions inside [0, 1]";
    } else {
        return true;
    }
    return false;
}

// ============================================================================
// Construction & Initialization
// ============================================================================

FLIPSimulator::FLIPSimulator(const SimConfig& cfg)
    : config(checkedConfig(cfg))
    , grid(config.gridNX, config.gridNY, config.cellSize())
    , pressureSolver(config.pressureIterations, config.density, config.pressureTolerance)
{
}

void FLIPSimulator::initialize() {
    initializeGrid();
    initializeParticles();
    std::cout << "[Init] Grid: " << config.gridNX << "x" << config.gridNY
              << ", dx=" << config.cellSize()
              << ", domain=" << config.domainWidth << "x" << config.domainHeight()
              << std::endl;
    std::cout << "[Init] Particles: " << particles.size()
              << ", FLIP alpha=" << config.flipAlpha
              << ", dt=" << config.dt
              << std::endl;
}

void FLIPSimulator::initializeGrid() {
    // Mark boundary cells as SOLID (one-cell-thick border)
    for (int j = 0; j < grid.NY; ++j) {
        for (int i = 0; i < grid.NX; ++i) {
            if (i == 0 || i == grid.NX - 1 || j == 0 || j == grid.NY - 1) {
                grid.cellType[grid.cellIndex(i, j)] = CellType::SOLID;
            } else {
                grid.cellType[grid.cellIndex(i, j)] = CellType::AIR;
            }
        }
    }
}

void FLIPSimulator::initializeParticles() {
    particles.clear();

    float dx = config.cellSize();
    float domH = config.domainHeight();

    // Dam-break region in world coordinates (offset by 1 cell to avoid solids)
    float xMin = config.damLeft * config.domainWidth + dx;   // stay inside solid border
    float xMax = config.damRight * config.domainWidth;
    float yMin = config.damBottom * domH + dx;
    float yMax = config.damTop * domH;

    // Clamp to be inside the solid boundary
    float solidMargin = dx + config.particleRadius;
    xMin = std::max(xMin, solidMargin);
    xMax = std::min(xMax, config.domainWidth - solidMargin);
    yMin = std::max(yMin, solidMargin);
    yMax = std::min(yMax, domH - solidMargin);

    // Place particles in a sub-grid pattern within each cell
    // particlesPerCell = 4 means 2x2 sub-grid
    int subN = static_cast<int>(std::sqrt(static_cast<float>(config.particlesPerCell)));
    if (subN < 1) subN = 1;
    float subDx = dx / static_cast<float>(subN);

    std::mt19937 rng(config.seed);
    std::uniform_real_distribution<float> jitter(-0.25f * subDx, 0.25f * subDx);

    for (float y = yMin; y < yMax; y += subDx) {
        for (float x = xMin; x < xMax; x += subDx) {
            float px = x + 0.5f * subDx + jitter(rng);
            float py = y + 0.5f * subDx + jitter(rng);

            // Ensure particle is inside domain (not in solid cells)
            px = std::max(solidMargin, std::min(px, config.domainWidth - solidMargin));
            py = std::max(solidMargin, std::min(py, domH - solidMargin));

            // Verify the cell is not SOLID
            int ci = static_cast<int>(px / dx);
            int cj = static_cast<int>(py / dx);
            ci = std::max(0, std::min(ci, grid.NX - 1));
            cj = std::max(0, std::min(cj, grid.NY - 1));
            if (grid.cellType[grid.cellIndex(ci, cj)] == CellType::SOLID) {
                continue;
            }

            Particle p;
            p.x = px;
            p.y = py;
            p.vx = 0.0f;
            p.vy = 0.0f;
            particles.push_back(p);
        }
    }

    std::cout << "[Init] Created " << particles.size() << " particles in dam-break region."
              << std::endl;
}

// ============================================================================
// Main simulation loop
// ============================================================================

void FLIPSimulator::step(float dt) {
    // STAGE 2: Particle to Grid
    particleToGrid();

    // Classify cells
    updateCellTypes();

    extrapolateVelocities();

    // Snapshot velocities for FLIP (before forces and pressure)
    savePreviousVelocities();

    // STAGE 3: External forces
    applyGravity(dt);

    // Enforce solid boundary conditions
    applyBoundaryConditions();

    // STAGE 4: Pressure projection
    solvePressure(dt);

    // Re-enforce boundary conditions after pressure correction
    applyBoundaryConditions();

    pressureSolver.computeDivergence(grid);
    divAfterProjection = pressureSolver.computeMaxDivergence(grid);

    // STAGE 5: Grid to Particle
    gridToParticle();

    // STAGE 6: Advection and collisions
    updateParticlePositions(dt);
    handleCollisions();
}

void FLIPSimulator::run() {
    if (config.enableOutput && config.outputInterval > 0) {
    prepareParticleOutput(config.outputDir);
    }

    for (int s = 0; s < config.numSteps; ++s) {
        currentStep = s;
        step(config.dt);

        // Output particle data at configured intervals
        if (config.enableOutput && config.outputInterval > 0 &&
            (s % config.outputInterval == 0 || s == config.numSteps - 1)) {
            writeParticleData(s);
        }

        // Progress report
        if ((s + 1) % 50 == 0 || s == 0 || s == config.numSteps - 1) {
            std::cout << "[Step " << std::setw(4) << (s + 1) << "/" << config.numSteps
                      << "] particles=" << particles.size()
                      << "  divBefore=" << std::scientific << std::setprecision(4) << divBeforeProjection
                      << "  divAfter=" << divAfterProjection
                      << std::fixed << std::setprecision(4)
                      << std::endl;
        }
    }
}

// ============================================================================
// STAGE 2: Particle to Grid (P2G)
// ============================================================================

void FLIPSimulator::particleToGrid() {
    float dx = grid.dx;
    float invDx = grid.invDx;
    int NX = grid.NX;
    int NY = grid.NY;

    // Clear velocity and weight accumulators
    std::fill(grid.u.begin(), grid.u.end(), 0.0f);
    std::fill(grid.uWeight.begin(), grid.uWeight.end(), 0.0f);
    std::fill(grid.v.begin(), grid.v.end(), 0.0f);
    std::fill(grid.vWeight.begin(), grid.vWeight.end(), 0.0f);

    for (const auto& p : particles) {
        // ---- Scatter vx onto u-faces ----
        // u-faces are at positions (i*dx, (j+0.5)*dx)
        // Convert particle position to u-face grid coordinates
        float ux = p.x * invDx;
        float uy = p.y * invDx - 0.5f;

        int i0 = static_cast<int>(std::floor(ux));
        int j0 = static_cast<int>(std::floor(uy));

        float sx = ux - static_cast<float>(i0);  // fractional part in x
        float sy = uy - static_cast<float>(j0);  // fractional part in y

        // Bilinear weights
        float w00 = (1.0f - sx) * (1.0f - sy);
        float w10 = sx * (1.0f - sy);
        float w01 = (1.0f - sx) * sy;
        float w11 = sx * sy;

        // Accumulate onto the four nearest u-faces
        // u has dimensions (NX+1) x NY, valid indices: i ∈ [0,NX], j ∈ [0,NY-1]
        for (int dj = 0; dj <= 1; ++dj) {
            for (int di = 0; di <= 1; ++di) {
                int ii = i0 + di;
                int jj = j0 + dj;
                if (ii >= 0 && ii <= NX && jj >= 0 && jj < NY) {
                    float w = (di == 0 ? (1.0f - sx) : sx) * (dj == 0 ? (1.0f - sy) : sy);
                    int idx = grid.uIndex(ii, jj);
                    grid.u[idx] += w * p.vx;
                    grid.uWeight[idx] += w;
                }
            }
        }

        // ---- Scatter vy onto v-faces ----
        // v-faces are at positions ((i+0.5)*dx, j*dx)
        float vx = p.x * invDx - 0.5f;
        float vy = p.y * invDx;

        i0 = static_cast<int>(std::floor(vx));
        j0 = static_cast<int>(std::floor(vy));

        sx = vx - static_cast<float>(i0);
        sy = vy - static_cast<float>(j0);

        // Accumulate onto the four nearest v-faces
        // v has dimensions NX x (NY+1), valid indices: i ∈ [0,NX-1], j ∈ [0,NY]
        for (int dj = 0; dj <= 1; ++dj) {
            for (int di = 0; di <= 1; ++di) {
                int ii = i0 + di;
                int jj = j0 + dj;
                if (ii >= 0 && ii < NX && jj >= 0 && jj <= NY) {
                    float w = (di == 0 ? (1.0f - sx) : sx) * (dj == 0 ? (1.0f - sy) : sy);
                    int idx = grid.vIndex(ii, jj);
                    grid.v[idx] += w * p.vy;
                    grid.vWeight[idx] += w;
                }
            }
        }
    }

    // Normalize by weights
    for (int idx = 0; idx < grid.uSize(); ++idx) {
        if (grid.uWeight[idx] > 0.0f) {
            grid.u[idx] /= grid.uWeight[idx];
        } else {
            grid.u[idx] = 0.0f;
        }
    }
    for (int idx = 0; idx < grid.vSize(); ++idx) {
        if (grid.vWeight[idx] > 0.0f) {
            grid.v[idx] /= grid.vWeight[idx];
        } else {
            grid.v[idx] = 0.0f;
        }
    }
}

// ============================================================================
// Cell type classification
// ============================================================================

void FLIPSimulator::updateCellTypes() {
    float dx = grid.dx;
    int NX = grid.NX;
    int NY = grid.NY;

    // Reset non-solid cells to AIR
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = grid.cellIndex(i, j);
            if (grid.cellType[idx] != CellType::SOLID) {
                grid.cellType[idx] = CellType::AIR;
            }
        }
    }

    // Mark cells containing particles as FLUID
    for (const auto& p : particles) {
        int ci = static_cast<int>(p.x / dx);
        int cj = static_cast<int>(p.y / dx);
        ci = std::max(0, std::min(ci, NX - 1));
        cj = std::max(0, std::min(cj, NY - 1));
        int idx = grid.cellIndex(ci, cj);
        if (grid.cellType[idx] != CellType::SOLID) {
            grid.cellType[idx] = CellType::FLUID;
        }
    }
}

void FLIPSimulator::extrapolateVelocities() {
    for (int layer = 0; layer < 3; ++layer) {
    for (int j = 0; j < grid.NY; ++j) {
        for (int i = 0; i <= grid.NX; ++i) {
                const int face = grid.uIndex(i, j);
                if (grid.uWeight[face] != 0.0f) {
                    continue;
                }

                const bool leftFluid = i > 0 &&
                    grid.cellType[grid.cellIndex(i - 1, j)] == CellType::FLUID;
                const bool rightFluid = i < grid.NX &&
                    grid.cellType[grid.cellIndex(i, j)] == CellType::FLUID;
                if (!leftFluid && !rightFluid) {
                    continue;
                }

                float sum = 0.0f;
                int count = 0;
                for (int offset : {-1, 1}) {
                    const int neighbor = i + offset;
                    if (neighbor >= 0 && neighbor <= grid.NX) {
                        const int neighborFace = grid.uIndex(neighbor, j);
                        if (grid.uWeight[neighborFace] > 0.0f) {
                            sum += grid.u[neighborFace];
                            ++count;
                        }
                    }
                }
                for (int offset : {-1, 1}) {
                    const int neighbor = j + offset;
                    if (neighbor >= 0 && neighbor < grid.NY) {
                        const int neighborFace = grid.uIndex(i, neighbor);
                        if (grid.uWeight[neighborFace] > 0.0f) {
                            sum += grid.u[neighborFace];
                            ++count;
                        }
                    }
                }
                if (count > 0) {
                    grid.u[face] = sum / static_cast<float>(count);
                    grid.uWeight[face] = -1.0f;
                }
            }
        }

    for (int j = 0; j <= grid.NY; ++j) {
        for (int i = 0; i < grid.NX; ++i) {
                const int face = grid.vIndex(i, j);
                if (grid.vWeight[face] != 0.0f) {
                    continue;
                }

                const bool belowFluid = j > 0 &&
                    grid.cellType[grid.cellIndex(i, j - 1)] == CellType::FLUID;
                const bool aboveFluid = j < grid.NY &&
                    grid.cellType[grid.cellIndex(i, j)] == CellType::FLUID;
                if (!belowFluid && !aboveFluid) {
                    continue;
                }

                float sum = 0.0f;
                int count = 0;
                for (int offset : {-1, 1}) {
                    const int neighbor = i + offset;
                    if (neighbor >= 0 && neighbor < grid.NX) {
                        const int neighborFace = grid.vIndex(neighbor, j);
                        if (grid.vWeight[neighborFace] > 0.0f) {
                            sum += grid.v[neighborFace];
                            ++count;
                        }
                    }
                }
                for (int offset : {-1, 1}) {
                    const int neighbor = j + offset;
                    if (neighbor >= 0 && neighbor <= grid.NY) {
                        const int neighborFace = grid.vIndex(i, neighbor);
                        if (grid.vWeight[neighborFace] > 0.0f) {
                            sum += grid.v[neighborFace];
                            ++count;
                        }
                    }
                }
                if (count > 0) {
                    grid.v[face] = sum / static_cast<float>(count);
                    grid.vWeight[face] = -1.0f;
                }
            }
        }
    }

    for (float& weight : grid.uWeight) {
        if (weight < 0.0f) weight = 1.0f;
    }
    for (float& weight : grid.vWeight) {
        if (weight < 0.0f) weight = 1.0f;
    }
}

// ============================================================================
// Save previous velocities (for FLIP)
// ============================================================================

void FLIPSimulator::savePreviousVelocities() {
    // Copy current grid velocities to "previous" buffers.
    // These are the post-P2G velocities BEFORE forces and pressure correction.
    // The FLIP update computes: delta = corrected - previous
    std::copy(grid.u.begin(), grid.u.end(), grid.uPrev.begin());
    std::copy(grid.v.begin(), grid.v.end(), grid.vPrev.begin());
}

// ============================================================================
// STAGE 3: External forces (gravity)
// ============================================================================

void FLIPSimulator::applyGravity(float dt) {
    int NX = grid.NX;
    int NY = grid.NY;

    // Add gravity to v-faces (vertical velocity).
    // Only apply to faces that are adjacent to at least one FLUID cell.
    // v[i,j] is the vertical face between cell (i,j-1) and (i,j).
    for (int j = 1; j < NY; ++j) {   // j=0 and j=NY are domain boundaries
        for (int i = 0; i < NX; ++i) {
            int cellBelow = grid.cellIndex(i, j - 1);
            int cellAbove = grid.cellIndex(i, j);
            // Apply gravity if at least one adjacent cell is FLUID
            if (grid.cellType[cellBelow] == CellType::FLUID ||
                grid.cellType[cellAbove] == CellType::FLUID) {
                grid.v[grid.vIndex(i, j)] += config.gravity * dt;
            }
        }
    }
}

// ============================================================================
// Boundary conditions
// ============================================================================

void FLIPSimulator::applyBoundaryConditions() {
    int NX = grid.NX;
    int NY = grid.NY;

    // Zero out velocity components at solid boundaries.
    // For u-faces: u[i,j] is between cell (i-1,j) and (i,j).
    //   If either neighbor is SOLID, the face velocity should not allow flow into the solid.
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i <= NX; ++i) {
            // Left or right neighbor is SOLID → zero normal velocity
            bool leftSolid  = (i > 0)  ? (grid.cellType[grid.cellIndex(i - 1, j)] == CellType::SOLID) : true;
            bool rightSolid = (i < NX) ? (grid.cellType[grid.cellIndex(i, j)]     == CellType::SOLID) : true;
            if (leftSolid || rightSolid) {
                grid.u[grid.uIndex(i, j)] = 0.0f;
            }
        }
    }

    // For v-faces: v[i,j] is between cell (i,j-1) and (i,j).
    for (int j = 0; j <= NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            bool belowSolid = (j > 0)  ? (grid.cellType[grid.cellIndex(i, j - 1)] == CellType::SOLID) : true;
            bool aboveSolid = (j < NY) ? (grid.cellType[grid.cellIndex(i, j)]     == CellType::SOLID) : true;
            if (belowSolid || aboveSolid) {
                grid.v[grid.vIndex(i, j)] = 0.0f;
            }
        }
    }

        // Tangencijalne brzine uz zidove (free-slip)
    for (int i = 1; i < NX; ++i) {
        grid.u[grid.uIndex(i, 0)] =
            grid.u[grid.uIndex(i, 1)];

        grid.u[grid.uIndex(i, NY - 1)] =
            grid.u[grid.uIndex(i, NY - 2)];
    }

    for (int j = 1; j < NY; ++j) {
        grid.v[grid.vIndex(0, j)] =
            grid.v[grid.vIndex(1, j)];

        grid.v[grid.vIndex(NX - 1, j)] =
            grid.v[grid.vIndex(NX - 2, j)];
    }
}

// ============================================================================
// STAGE 4: Pressure projection
// ===========================================================================

void FLIPSimulator::solvePressure(float dt) {
    // Compute divergence before projection (for reporting)
    pressureSolver.computeDivergence(grid);
    divBeforeProjection = pressureSolver.computeMaxDivergence(grid);

    // Solve pressure and correct velocities
    pressureSolver.solve(grid, dt);

    // Compute divergence after projection (for reporting)
    pressureSolver.computeDivergence(grid);
    divAfterProjection = pressureSolver.computeMaxDivergence(grid);
}

// ============================================================================
// STAGE 5: Grid to Particle (G2P)
// ============================================================================

void FLIPSimulator::gridToParticle() {
    float alpha = config.flipAlpha;

    for (auto& p : particles) {
        // PIC: interpolate corrected grid velocity directly
        float picVx = interpolateU(p.x, p.y, grid.u);
        float picVy = interpolateV(p.x, p.y, grid.v);

        // FLIP: particle velocity + grid velocity change
        float prevVx = interpolateU(p.x, p.y, grid.uPrev);
        float prevVy = interpolateV(p.x, p.y, grid.vPrev);
        float flipVx = p.vx + (picVx - prevVx);
        float flipVy = p.vy + (picVy - prevVy);

        // Blend PIC and FLIP
        p.vx = (1.0f - alpha) * picVx + alpha * flipVx;
        p.vy = (1.0f - alpha) * picVy + alpha * flipVy;
    }
}

// ============================================================================
// Bilinear interpolation of staggered velocities
// ============================================================================

float FLIPSimulator::interpolateU(float px, float py, const std::vector<float>& uField) const {
    // u-faces are at positions (i*dx, (j+0.5)*dx)
    // Convert world position to u-face grid coordinates
    float dx = grid.dx;
    float invDx = grid.invDx;

    float gx = px * invDx;          // continuous grid x-index
    float gy = py * invDx - 0.5f;   // offset because u is at cell-center height

    int i0 = static_cast<int>(std::floor(gx));
    int j0 = static_cast<int>(std::floor(gy));

    float sx = gx - static_cast<float>(i0);
    float sy = gy - static_cast<float>(j0);

    // Clamp indices to valid range: u has (NX+1) columns and NY rows
    auto clampI = [&](int i) { return std::max(0, std::min(i, grid.NX)); };
    auto clampJ = [&](int j) { return std::max(0, std::min(j, grid.NY - 1)); };

    int i0c = clampI(i0);
    int i1c = clampI(i0 + 1);
    int j0c = clampJ(j0);
    int j1c = clampJ(j0 + 1);

    float v00 = uField[grid.uIndex(i0c, j0c)];
    float v10 = uField[grid.uIndex(i1c, j0c)];
    float v01 = uField[grid.uIndex(i0c, j1c)];
    float v11 = uField[grid.uIndex(i1c, j1c)];

    return (1.0f - sx) * (1.0f - sy) * v00
         +        sx  * (1.0f - sy) * v10
         + (1.0f - sx) *        sy  * v01
         +        sx  *        sy  * v11;
}

float FLIPSimulator::interpolateV(float px, float py, const std::vector<float>& vField) const {
    // v-faces are at positions ((i+0.5)*dx, j*dx)
    float dx = grid.dx;
    float invDx = grid.invDx;

    float gx = px * invDx - 0.5f;   // offset because v is at cell-center width
    float gy = py * invDx;          // continuous grid y-index

    int i0 = static_cast<int>(std::floor(gx));
    int j0 = static_cast<int>(std::floor(gy));

    float sx = gx - static_cast<float>(i0);
    float sy = gy - static_cast<float>(j0);

    // Clamp indices: v has NX columns and (NY+1) rows
    auto clampI = [&](int i) { return std::max(0, std::min(i, grid.NX - 1)); };
    auto clampJ = [&](int j) { return std::max(0, std::min(j, grid.NY)); };

    int i0c = clampI(i0);
    int i1c = clampI(i0 + 1);
    int j0c = clampJ(j0);
    int j1c = clampJ(j0 + 1);

    float v00 = vField[grid.vIndex(i0c, j0c)];
    float v10 = vField[grid.vIndex(i1c, j0c)];
    float v01 = vField[grid.vIndex(i0c, j1c)];
    float v11 = vField[grid.vIndex(i1c, j1c)];

    return (1.0f - sx) * (1.0f - sy) * v00
         +        sx  * (1.0f - sy) * v10
         + (1.0f - sx) *        sy  * v01
         +        sx  *        sy  * v11;
}

// ============================================================================
// STAGE 6: Particle advection (RK2 midpoint)
// ============================================================================

void FLIPSimulator::updateParticlePositions(float dt) {
    // Use RK2 (midpoint) advection with the corrected grid velocity field.
    // This is more accurate than forward Euler and avoids excessive numerical diffusion.
    //
    // Step 1: Compute midpoint position
    //   x_mid = x + 0.5 * dt * v(x)
    // Step 2: Evaluate velocity at midpoint
    //   v_mid = interpolate(grid, x_mid)
    // Step 3: Full step
    //   x_new = x + dt * v_mid

    float minX = config.cellSize() + config.particleRadius;
    float maxX = config.domainWidth - config.cellSize() - config.particleRadius;
    float minY = config.cellSize() + config.particleRadius;
    float maxY = config.domainHeight() - config.cellSize() - config.particleRadius;

    for (auto& p : particles) {
        float maxSpeed = std::max(std::fabs(p.vx), std::fabs(p.vy));
        int substeps = std::max(1, static_cast<int>(std::ceil(
            maxSpeed * dt / (0.5f * grid.dx))));
        float subDt = dt / static_cast<float>(substeps);

        for (int substep = 0; substep < substeps; ++substep) {
            float vx0 = interpolateU(p.x, p.y, grid.u);
            float vy0 = interpolateV(p.x, p.y, grid.v);
            float mx = std::max(minX, std::min(p.x + 0.5f * subDt * vx0, maxX));
            float my = std::max(minY, std::min(p.y + 0.5f * subDt * vy0, maxY));
            float vxm = interpolateU(mx, my, grid.u);
            float vym = interpolateV(mx, my, grid.v);

            p.x += subDt * vxm;
            p.y += subDt * vym;

            if (p.x < minX) {
                p.x = minX;
                p.vx = 0.0f;
            } else if (p.x > maxX) {
                p.x = maxX;
                p.vx = 0.0f;
            }
            if (p.y < minY) {
                p.y = minY;
                p.vy = 0.0f;
            } else if (p.y > maxY) {
                p.y = maxY;
                p.vy = 0.0f;
            }
        }
    }
}

// ============================================================================
// Collision handling
// ============================================================================

void FLIPSimulator::handleCollisions() {
    float dx = grid.dx;
    float domW = config.domainWidth;
    float domH = config.domainHeight();

    // Boundary clearance: particle should stay at least one solid cell + radius away
    float minX = dx + config.particleRadius;
    float maxX = domW - dx - config.particleRadius;
    float minY = dx + config.particleRadius;
    float maxY = domH - dx - config.particleRadius;

    for (auto& p : particles) {
        // Left wall
        if (p.x < minX) {
            p.x = minX;
            if (p.vx < 0.0f) p.vx = 0.0f;
        }
        // Right wall
        if (p.x > maxX) {
            p.x = maxX;
            if (p.vx > 0.0f) p.vx = 0.0f;
        }
        // Bottom wall
        if (p.y < minY) {
            p.y = minY;
            if (p.vy < 0.0f) p.vy = 0.0f;
        }
        // Top wall
        if (p.y > maxY) {
            p.y = maxY;
            if (p.vy > 0.0f) p.vy = 0.0f;
        }
    }
}

// ============================================================================
// File output
// ============================================================================

void FLIPSimulator::writeParticleData(int stepNumber) const {
    std::ostringstream filename;
    filename << config.outputDir << "/particles_"
             << std::setw(5) << std::setfill('0') << stepNumber << ".csv";

    std::ofstream file(filename.str());
    if (!file.is_open()) {
        std::cerr << "[Warning] Could not open " << filename.str() << " for writing." << std::endl;
        return;
    }

    file << "x,y,vx,vy\n";
    file << std::fixed << std::setprecision(6);
    for (const auto& p : particles) {
        file << p.x << "," << p.y << "," << p.vx << "," << p.vy << "\n";
    }
}
