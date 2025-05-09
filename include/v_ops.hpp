#pragma once

#include <vector>
#include <stdexcept>

/**
 * Calculates the squared euclidean distance of N Dimensions, between two points.
 * @param p1
 * @param p2
 * @returns Squared euclidean distance.
 */
inline double sqrDist(const float *v1, const float *v2, uint32_t dim)
{
    double sqr_dist = 0;
    // #pragma omp simd //?
    //#pragma omp parallel for reduction(+:sqr_dist)
    for (int i = 0; i < dim; i++)
    {
        sqr_dist += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    
    return sqr_dist;
};

inline void addInto(double* v1, const float* v2, uint32_t dim)
{
    //#pragma omp simd
    for (int i = 0; i < dim; i++)
    {
        v1[i] += v2[i];
    }
}

inline void addInto(double* v1, const double* v2, uint32_t dim)
{
    //#pragma omp simd
    for (int i = 0; i < dim; i++)
    {
        v1[i] += v2[i];
    }
}
