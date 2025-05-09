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
    std::vector<double> add_rem_points_sums;
    std::vector<uint32_t> add_rem_points_counts;

    inline CentroidDiffs(int n, uint32_t dim) : add_rem_points_sums(n*dim*2, 0), add_rem_points_counts(n*2,0), n(n), dim(dim) {}

    inline auto getDim() const { return dim; }
    inline CentroidDiffRef operator[](int i) { return CentroidDiffRef{&add_rem_points_sums[i*dim], &add_rem_points_sums[n*dim+i*dim], &add_rem_points_counts[i], &add_rem_points_counts[n+i], dim}; }
};
