#ifndef PRESSURE_SOLVER_HPP
#define PRESSURE_SOLVER_HPP

/**
 * @file PressureSolver.hpp
 * @brief Pressure projection via Jacobi iteration on the MAC grid.
 *
 * Solves the discrete pressure Poisson equation:
 *
 *     ∑ A_ij * p_neighbor - A_diag * p_center = -divergence / dt
 *
 * where A coefficients account for FLUID/AIR/SOLID neighbor types.
 *
 * Uses double-buffered Jacobi iteration (reads from pressureOld,
 * writes to pressureNew, then swaps). This is trivially parallelizable
 * with OpenMP and vectorizable with SIMD in future versions.
 */

#include "Grid.hpp"

/**
 * @class PressureSolver
 * @brief Jacobi-based pressure solver for the incompressible Navier-Stokes equations.
 */
class PressureSolver {
public:
    int maxIterations;    ///< Maximum number of Jacobi iterations
    float density;        ///< Fluid density (kg/m³)
    float tolerance;      ///< Relative pressure-update tolerance
    mutable int lastIterations;
    mutable float lastResidual;
    mutable float lastPressureDelta;
    mutable float initialResidual;
    mutable float residualRatio;

    PressureSolver();
    PressureSolver(int maxIter, float rho, float convergenceTolerance = 1.0e-4f);

    /**
     * @brief Compute the velocity divergence field for all FLUID cells.
     * @param grid The MAC grid with current velocity fields.
     *
     * Divergence for cell (i,j):
     *   div = (u[i+1,j] - u[i,j] + v[i,j+1] - v[i,j]) / dx
     */
    void computeDivergence(Grid& grid) const;

    /**
     * @brief Solve the pressure Poisson equation using Jacobi iteration.
     * @param grid The MAC grid.
     * @param dt   Timestep.
     *
     * After solving, corrects face velocities using the pressure gradient:
     *   u[i,j] -= dt / (density * dx) * (p[i,j] - p[i-1,j])
     *   v[i,j] -= dt / (density * dx) * (p[i,j] - p[i,j-1])
     */
    void solve(Grid& grid, float dt) const;

    /**
     * @brief Apply pressure gradient to correct velocities.
     * @param grid The MAC grid with solved pressure field.
     * @param dt   Timestep.
     */
    void applyPressureGradient(Grid& grid, float dt) const;

    /**
     * @brief Compute and return the maximum absolute divergence across all FLUID cells.
     * @param grid The MAC grid.
     * @return Maximum |divergence| value.
     */
    float computeMaxDivergence(const Grid& grid) const;

    /**
     * @brief Compute and return the average absolute divergence across all FLUID cells.
     */
    float computeAvgDivergence(const Grid& grid) const;
};

#endif // PRESSURE_SOLVER_HPP
