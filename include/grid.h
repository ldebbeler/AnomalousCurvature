//File: grid.h
#ifndef GRID_H
#define GRID_H
#include <vector>

struct grid{
    std::vector<double> m_x;
    grid(double init, double fin, std::size_t L);
    std::vector<double> logCreate(double init, double fin, std::size_t L);
    std::vector<double> linCreate(double fin, std::size_t L);
};

#endif //GRID_H
