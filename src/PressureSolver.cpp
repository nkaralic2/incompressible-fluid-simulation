/**
 * @file PressureSolver.cpp
 * @brief Jacobi-based pressure solver for the incompressible Navier-Stokes equations.
 *
 * Solves the pressure Poisson equation on the MAC grid:
 *
 *   ∇²p = (ρ / dt) * ∇·u
 *
 * In discrete form for cell (i,j):
 *
 *   Σ_neighbors (p_neighbor - p_center) = -(dx² * ρ / dt) * div(i,j)
 *
 * where the sum is over non-solid neighbors, and:
 *   - AIR neighbors contribute p = 0 (atmospheric pressure).
 *   - FLUID neighbors contribute their unknown pressure.
 *   - SOLID neighbors are excluded from the sum (no flow through walls).
 *
 * The Jacobi update for cell (i,j):
 *
 *   p_new(i,j) = (Σ p_old(neighbor) + scale * div(i,j)) / numFluidNeighbors
 *
 * where scale = -(dx² * ρ / dt) and the sum includes only FLUID/AIR neighbors.
 *
 * After convergence, velocities are corrected:
 *   u(i,j) -= (dt / (ρ * dx)) * (p(i,j) - p(i-1,j))
 *   v(i,j) -= (dt / (ρ * dx)) * (p(i,j) - p(i,j-1))
 */

#include "PressureSolver.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

PressureSolver::PressureSolver()
        : maxIterations(100), density(1000.0f), tolerance(1.0e-4f),
                lastIterations(0), lastResidual(0.0f), lastPressureDelta(0.0f),
                initialResidual(0.0f), residualRatio(0.0f)
{
}

PressureSolver::PressureSolver(int maxIter, float rho, float convergenceTolerance)
        : maxIterations(maxIter), density(rho), tolerance(convergenceTolerance),
                lastIterations(0), lastResidual(0.0f), lastPressureDelta(0.0f),
                initialResidual(0.0f), residualRatio(0.0f)
{
}

// ============================================================================
// Divergence computation
// ============================================================================

void PressureSolver::computeDivergence(Grid& grid) const {
    int NX = grid.NX;
    int NY = grid.NY;
    float invDx = grid.invDx;

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = grid.cellIndex(i, j);
            if (grid.cellType[idx] == CellType::FLUID) {
                // Divergence = (u_right - u_left + v_top - v_bottom) / dx
                float uRight = grid.u[grid.uIndex(i + 1, j)];
                float uLeft  = grid.u[grid.uIndex(i, j)];
                float vTop   = grid.v[grid.vIndex(i, j + 1)];
                float vBot   = grid.v[grid.vIndex(i, j)];
                grid.divergence[idx] = (uRight - uLeft + vTop - vBot) * invDx;
            } else {
                grid.divergence[idx] = 0.0f;
            }
        }
    }
}

// ============================================================================
// Jacobi pressure solve
// ============================================================================

void PressureSolver::solve(Grid& grid, float dt) const {
    int NX = grid.NX;
    int NY = grid.NY;
    float dx = grid.dx;

    // Compute divergence of the current velocity field
    computeDivergence(grid);

    // Scale factor for the RHS of the Poisson equation
    // From: ∇²p = (ρ/dt) * ∇·u
    // Discrete: Σ(p_neighbor - p_center) = (dx² * ρ / dt) * div
    // Rearranged for Jacobi: p_center = (Σ p_neighbor - scale * div) / diag
    // where scale = dx² * ρ / dt
    float scale = dx * dx * density / dt;

    // Initialize pressure buffers
    std::fill(grid.pressureOld.begin(), grid.pressureOld.end(), 0.0f);
    std::fill(grid.pressureNew.begin(), grid.pressureNew.end(), 0.0f);

    lastIterations = maxIterations;
    lastResidual = 0.0f;
    lastPressureDelta = 0.0f;
    initialResidual = 0.0f;
    residualRatio = 0.0f;

    for (int j = 1; j < NY - 1; ++j) {
        for (int i = 1; i < NX - 1; ++i) {
            int idx = grid.cellIndex(i, j);
            if (grid.cellType[idx] == CellType::FLUID) {
                initialResidual = std::max(initialResidual,
                    std::fabs(scale * grid.divergence[idx]));
            }
        }
    }

    // Jacobi iteration with an adaptive stopping criterion based on the
    // relative pressure update.
    for (int iter = 0; iter < maxIterations; ++iter) {
        for (int j = 1; j < NY - 1; ++j) {
            for (int i = 1; i < NX - 1; ++i) {
                int idx = grid.cellIndex(i, j);
                if (grid.cellType[idx] != CellType::FLUID) {
                    grid.pressureNew[idx] = 0.0f;
                    continue;
                }

                // Count non-solid neighbors and accumulate their pressures
                float sumP = 0.0f;
                int numNeighbors = 0;  // number of non-solid neighbors (diagonal coefficient)

                // Left neighbor (i-1, j)
                {
                    int ni = i - 1;
                    CellType ct = grid.cellType[grid.cellIndex(ni, j)];
                    if (ct != CellType::SOLID) {
                        numNeighbors++;
                        if (ct == CellType::FLUID) {
                            sumP += grid.pressureOld[grid.cellIndex(ni, j)];
                        }
                        // AIR: pressure = 0, adds nothing
                    }
                }

                // Right neighbor (i+1, j)
                {
                    int ni = i + 1;
                    CellType ct = grid.cellType[grid.cellIndex(ni, j)];
                    if (ct != CellType::SOLID) {
                        numNeighbors++;
                        if (ct == CellType::FLUID) {
                            sumP += grid.pressureOld[grid.cellIndex(ni, j)];
                        }
                    }
                }

                // Bottom neighbor (i, j-1)
                {
                    int nj = j - 1;
                    CellType ct = grid.cellType[grid.cellIndex(i, nj)];
                    if (ct != CellType::SOLID) {
                        numNeighbors++;
                        if (ct == CellType::FLUID) {
                            sumP += grid.pressureOld[grid.cellIndex(i, nj)];
                        }
                    }
                }

                // Top neighbor (i, j+1)
                {
                    int nj = j + 1;
                    CellType ct = grid.cellType[grid.cellIndex(i, nj)];
                    if (ct != CellType::SOLID) {
                        numNeighbors++;
                        if (ct == CellType::FLUID) {
                            sumP += grid.pressureOld[grid.cellIndex(i, nj)];
                        }
                    }
                }

                if (numNeighbors > 0) {
                    // Jacobi update:
                    // p = (sumP - scale * div) / numNeighbors
                    grid.pressureNew[idx] = (sumP - scale * grid.divergence[idx]) /
                                            static_cast<float>(numNeighbors);
                } else {
                    grid.pressureNew[idx] = 0.0f;
                }
            }
        }

        float maxDelta = 0.0f;
        float maxPressure = 0.0f;
        float maxPoissonResidual = 0.0f;
        float maxRhs = 0.0f;
        for (int j = 1; j < NY - 1; ++j) {
            for (int i = 1; i < NX - 1; ++i) {
                int idx = grid.cellIndex(i, j);
                if (grid.cellType[idx] == CellType::FLUID) {
                    maxDelta = std::max(maxDelta,
                        std::fabs(grid.pressureNew[idx] - grid.pressureOld[idx]));
                    maxPressure = std::max(maxPressure, std::fabs(grid.pressureNew[idx]));

                    float neighborPressure = 0.0f;
                    int nonSolidNeighbors = 0;
                    const int neighborI[] = {i - 1, i + 1, i, i};
                    const int neighborJ[] = {j, j, j - 1, j + 1};
                    for (int neighbor = 0; neighbor < 4; ++neighbor) {
                        int neighborIdx = grid.cellIndex(neighborI[neighbor], neighborJ[neighbor]);
                        CellType neighborType = grid.cellType[neighborIdx];
                        if (neighborType != CellType::SOLID) {
                            ++nonSolidNeighbors;
                            if (neighborType == CellType::FLUID) {
                                neighborPressure += grid.pressureNew[neighborIdx];
                            }
                        }
                    }
                    float rhs = scale * grid.divergence[idx];
                    float residual = neighborPressure -
                        static_cast<float>(nonSolidNeighbors) * grid.pressureNew[idx] - rhs;
                    maxPoissonResidual = std::max(maxPoissonResidual, std::fabs(residual));
                    maxRhs = std::max(maxRhs, std::fabs(rhs));
                }
            }
        }

        std::swap(grid.pressureOld, grid.pressureNew);
        lastPressureDelta = maxDelta;
        lastResidual = maxPoissonResidual;
        residualRatio = initialResidual > 0.0f
            ? lastResidual / initialResidual : 0.0f;
        lastIterations = iter + 1;
        if (maxDelta <= tolerance * std::max(1.0f, maxPressure) &&
            maxPoissonResidual <= tolerance * std::max(1.0f, maxRhs)) {
            break;
        }
    }

    // Final result is in pressureOld after the last swap
    std::copy(grid.pressureOld.begin(), grid.pressureOld.end(), grid.pressure.begin());

    // Apply pressure gradient to correct velocities
    applyPressureGradient(grid, dt);
}

// ============================================================================
// Pressure gradient velocity correction
// ============================================================================

void PressureSolver::applyPressureGradient(Grid& grid, float dt) const {
    int NX = grid.NX;
    int NY = grid.NY;
    float dx = grid.dx;

    // Correction scale: dt / (ρ * dx)
    float gradScale = dt / (density * dx);

    // Correct u-velocities at interior vertical faces
    // u[i,j] is between cell (i-1,j) and (i,j)
    for (int j = 1; j < NY - 1; ++j) {
        for (int i = 1; i < NX; ++i) {
            int leftCell  = grid.cellIndex(i - 1, j);
            int rightCell = grid.cellIndex(i, j);

            // Only correct if at least one adjacent cell is FLUID
            // and neither is SOLID (boundary conditions handle SOLID faces)
            if (grid.cellType[leftCell] == CellType::SOLID ||
                grid.cellType[rightCell] == CellType::SOLID) {
                continue;
            }
            if (grid.cellType[leftCell] == CellType::FLUID ||
                grid.cellType[rightCell] == CellType::FLUID) {
                float pRight = grid.pressure[rightCell];
                float pLeft  = grid.pressure[leftCell];
                grid.u[grid.uIndex(i, j)] -= gradScale * (pRight - pLeft);
            }
        }
    }

    // Correct v-velocities at interior horizontal faces
    // v[i,j] is between cell (i,j-1) and (i,j)
    for (int j = 1; j < NY; ++j) {
        for (int i = 1; i < NX - 1; ++i) {
            int belowCell = grid.cellIndex(i, j - 1);
            int aboveCell = grid.cellIndex(i, j);

            if (grid.cellType[belowCell] == CellType::SOLID ||
                grid.cellType[aboveCell] == CellType::SOLID) {
                continue;
            }
            if (grid.cellType[belowCell] == CellType::FLUID ||
                grid.cellType[aboveCell] == CellType::FLUID) {
                float pAbove = grid.pressure[aboveCell];
                float pBelow = grid.pressure[belowCell];
                grid.v[grid.vIndex(i, j)] -= gradScale * (pAbove - pBelow);
            }
        }
    }
}

// ============================================================================
// Divergence metrics
// ============================================================================

float PressureSolver::computeMaxDivergence(const Grid& grid) const {
    float maxDiv = 0.0f;
    int NX = grid.NX;
    int NY = grid.NY;

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = grid.cellIndex(i, j);
            if (grid.cellType[idx] == CellType::FLUID) {
                float d = std::fabs(grid.divergence[idx]);
                if (d > maxDiv) maxDiv = d;
            }
        }
    }
    return maxDiv;
}

float PressureSolver::computeAvgDivergence(const Grid& grid) const {
    float sumDiv = 0.0f;
    int count = 0;
    int NX = grid.NX;
    int NY = grid.NY;

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = grid.cellIndex(i, j);
            if (grid.cellType[idx] == CellType::FLUID) {
                sumDiv += std::fabs(grid.divergence[idx]);
                count++;
            }
        }
    }
    return (count > 0) ? (sumDiv / static_cast<float>(count)) : 0.0f;
}
