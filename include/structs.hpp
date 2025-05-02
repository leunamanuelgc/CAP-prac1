#pragma once

#include "Point.hpp"

#include <cstdint>
#include <vector>


struct PointData
{
    uint32_t n_points;
    uint32_t n_dim;
    std::vector<Point> points;
};

struct CentroidDiff
{
    std::vector<double> add_points_sum;
    int add_points_count;
    std::vector<double> rem_points_sum;
    int rem_points_count;

    //inline CentroidDiff() {}
    inline CentroidDiff(uint32_t dim) : add_points_count(0), rem_points_count(0)
    {
        add_points_sum = std::vector<double>(dim, 0);
        rem_points_sum = std::vector<double>(dim, 0);
    }
};