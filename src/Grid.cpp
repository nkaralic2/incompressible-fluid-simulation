/**
 * @file Grid.cpp
 * @brief Implementation of the MAC staggered grid.
 */

#include "Grid.hpp"
#include <algorithm>

Grid::Grid() : NX(0), NY(0), dx(1.0f), invDx(1.0f) {}

Grid::Grid(int nx, int ny, float cellSize)
    : NX(nx), NY(ny), dx(cellSize), invDx(1.0f / cellSize)
{
    allocate();
}

void Grid::allocate() {
    int uSz = uSize();
    int vSz = vSize();
    int cSz = cellSize();

    u.resize(uSz, 0.0f);
    uPrev.resize(uSz, 0.0f);
    uWeight.resize(uSz, 0.0f);

    v.resize(vSz, 0.0f);
    vPrev.resize(vSz, 0.0f);
    vWeight.resize(vSz, 0.0f);

    pressure.resize(cSz, 0.0f);
    pressureOld.resize(cSz, 0.0f);
    pressureNew.resize(cSz, 0.0f);
    divergence.resize(cSz, 0.0f);

    cellType.resize(cSz, CellType::AIR);
}

void Grid::clear() {
    std::fill(u.begin(), u.end(), 0.0f);
    std::fill(uPrev.begin(), uPrev.end(), 0.0f);
    std::fill(uWeight.begin(), uWeight.end(), 0.0f);

    std::fill(v.begin(), v.end(), 0.0f);
    std::fill(vPrev.begin(), vPrev.end(), 0.0f);
    std::fill(vWeight.begin(), vWeight.end(), 0.0f);

    std::fill(pressure.begin(), pressure.end(), 0.0f);
    std::fill(pressureOld.begin(), pressureOld.end(), 0.0f);
    std::fill(pressureNew.begin(), pressureNew.end(), 0.0f);
    std::fill(divergence.begin(), divergence.end(), 0.0f);

    // Reset non-solid cells to AIR; keep SOLID cells
    for (auto& ct : cellType) {
        if (ct != CellType::SOLID) {
            ct = CellType::AIR;
        }
    }
}
