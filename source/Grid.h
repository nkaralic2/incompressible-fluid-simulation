
#ifndef GRID_H
#define GRID_H

#include <vector>

class Grid {
public:
    int nx, ny;
    double cellSize;

    // Horizontalne brzine (U)
    std::vector<double> u;

    // Vertikalne brzine (V)
    std::vector<double> v;

    // Pritisak u centrima celija
    std::vector<double> pressure;

    // Divergencija
    std::vector<double> divergence;

    // Oznake celija: fluid, zrak, prepreka
    std::vector<int> cellType;

    Grid(int nx, int ny, double cellSize);

    int index(int i, int j) const;

    int uIndex(int i, int j) const;

    int vIndex(int i, int j) const;

    void clear();
};

#endif
