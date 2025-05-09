#pragma once

#include "CentroidDiffRef.hpp"

#include <vector>
#include <cstdint>

struct CentroidDiffs
{
private:
    int n;
    uint32_t dim;
public:
    std::vector<double> add_points_sums;
    std::vector<double> rem_points_sums;
    std::vector<uint32_t> add_points_counts;
    std::vector<uint32_t> rem_points_counts;

    inline CentroidDiffRef operator[](int i) { return CentroidDiffRef{&add_points_sums[i*dim], &rem_points_sums[i*dim], &add_points_counts[i], &rem_points_counts[i], dim}; }

    inline CentroidDiffs(int n, uint32_t dim) : add_points_sums(n*dim, 0), rem_points_sums(n*dim, 0), add_points_counts(0), rem_points_counts(0), n(n), dim(dim) {}
};