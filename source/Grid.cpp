
#include "Grid.h"
#include <algorithm>

Grid::Grid(int nx, int ny, double cellSize)
    : nx(nx),
    ny(ny),
    cellSize(cellSize),
    u((nx + 1)* ny, 0.0),
    v(nx* (ny + 1), 0.0),
    pressure(nx* ny, 0.0),
    divergence(nx* ny, 0.0),
    cellType(nx* ny, 0)
{
}

int Grid::index(int i, int j) const {
    return i + j * nx;
}

int Grid::uIndex(int i, int j) const {
    return i + j * (nx + 1);
}

int Grid::vIndex(int i, int j) const {
    return i + j * nx;
}

void Grid::clear() {
    std::fill(u.begin(), u.end(), 0.0);
    std::fill(v.begin(), v.end(), 0.0);
    std::fill(pressure.begin(), pressure.end(), 0.0);
    std::fill(divergence.begin(), divergence.end(), 0.0);
}
