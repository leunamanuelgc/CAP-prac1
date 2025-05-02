#pragma once

#include <vector>
#include <stdexcept>

/**
 * Calculates the squared euclidean distance of N Dimensions, between two points.
 * @param p1
 * @param p2
 * @returns Squared euclidean distance.
 */
inline double sqrDist(const std::vector<float>& v1, const std::vector<float>& v2) {
    auto n_dim = v1.size();
    if (n_dim != v2.size())
        throw std::invalid_argument("Can't measure distance between points with different dimensions [P1-%d <=> P2-%d]");

    double sqr_dist = 0;
    #pragma omp parallel for default(shared) reduction(+:sqr_dist)
    for (int i = 0; i < n_dim; i++)
    {
        sqr_dist += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    
    return sqr_dist;
};
