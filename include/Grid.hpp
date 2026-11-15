#ifndef GRID_HPP
#define GRID_HPP

/**
 * @file Grid.hpp
 * @brief MAC (Marker-and-Cell) staggered grid for incompressible fluid simulation.
 *
 * Layout on a MAC grid:
 *   - Pressure (p) is stored at cell centers.
 *   - Horizontal velocity (u) is stored at vertical cell faces (left/right edges).
 *   - Vertical velocity (v) is stored at horizontal cell faces (bottom/top edges).
 *
 * For an NX x NY grid of cells:
 *   - u has dimensions (NX+1) x NY  — one extra column of vertical faces.
 *   - v has dimensions NX x (NY+1)  — one extra row of horizontal faces.
 *   - p has dimensions NX x NY.
 *
 * Indexing convention (row-major):
 *   - u[j * (NX+1) + i]  for face at column i, row j
 *   - v[j * NX + i]       for face at column i, row j
 *   - p[j * NX + i]       for cell at column i, row j
 *
 * Cell types:
 *   - SOLID: boundary cells (the outermost ring of cells).
 *   - FLUID: cells containing at least one particle.
 *   - AIR:   empty cells adjacent to fluid.
 */

#include <vector>
#include <cstdint>

/// Cell classification for the MAC grid.
enum class CellType : uint8_t {
    AIR   = 0,
    FLUID = 1,
    SOLID = 2
};

/**
 * @class Grid
 * @brief Stores all Eulerian field data on a 2D MAC staggered grid.
 *
 * All arrays use contiguous std::vector<float> storage for cache-friendliness
 * and future SIMD / OpenMP compatibility.
 */
class Grid {
public:
    int NX;        ///< Number of cells in X direction
    int NY;        ///< Number of cells in Y direction
    float dx;      ///< Cell size (uniform, square cells)
    float invDx;   ///< 1.0f / dx, cached for performance

    // --- Staggered velocity fields ---
    // u: horizontal velocity at vertical faces, dimensions (NX+1) x NY
    std::vector<float> u;
    std::vector<float> uPrev;     ///< Snapshot before forces/pressure for FLIP
    std::vector<float> uWeight;   ///< Accumulation weights for P2G

    // v: vertical velocity at horizontal faces, dimensions NX x (NY+1)
    std::vector<float> v;
    std::vector<float> vPrev;     ///< Snapshot before forces/pressure for FLIP
    std::vector<float> vWeight;   ///< Accumulation weights for P2G

    // --- Pressure solve ---
    std::vector<float> pressure;
    std::vector<float> pressureOld; ///< Jacobi double-buffer
    std::vector<float> pressureNew; ///< Jacobi double-buffer
    std::vector<float> divergence;  ///< Velocity divergence at cell centers

    // --- Cell classification ---
    std::vector<CellType> cellType;

    /// Construct a grid with the given cell resolution and cell size.
    Grid();
    Grid(int nx, int ny, float cellSize);

    /// Reset all velocity, pressure, divergence, and weight arrays to zero.
    void clear();

    /// (Re-)allocate all arrays to match current NX, NY.
    void allocate();

    // --- Index helpers (inline for performance) ---

    /// Linear index into u-array:  column i ∈ [0, NX], row j ∈ [0, NY-1]
    inline int uIndex(int i, int j) const { return j * (NX + 1) + i; }

    /// Linear index into v-array:  column i ∈ [0, NX-1], row j ∈ [0, NY]
    inline int vIndex(int i, int j) const { return j * NX + i; }

    /// Linear index into cell-centered arrays (pressure, divergence, cellType)
    inline int cellIndex(int i, int j) const { return j * NX + i; }

    /// Number of u-face entries
    inline int uSize() const { return (NX + 1) * NY; }

    /// Number of v-face entries
    inline int vSize() const { return NX * (NY + 1); }

    /// Number of cell-center entries
    inline int cellSize() const { return NX * NY; }
};

#endif // GRID_HPP
