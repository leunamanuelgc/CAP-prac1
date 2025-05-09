#pragma once

#include "PointRef.hpp"

#define NONE_CLUSTER -1

#include <cstdint>
#include <vector>
#include <mpi.h>
#include <memory>

// Struct of Arrays approach
struct PointData
{
private:
    uint32_t n_total_points;
    uint32_t n_dim;
public:
    std::vector<float> coords;
    std::vector<int> cluster_ids;

    PointData() = default;
    inline PointData(uint32_t n_points, uint32_t n_total_points, uint32_t n_dim) :
        coords(n_points * n_dim, 0),
        cluster_ids(n_points, NONE_CLUSTER),
        n_total_points(n_total_points),
        n_dim(n_dim) {}

    inline float* getPointCoords(uint32_t idx) { return &coords[idx * n_dim]; }
    inline int* getPointClusterID(uint32_t idx) { return &cluster_ids[idx]; }
    inline uint32_t getDim() const { return n_dim; }
    inline uint32_t getNPoints() const { return cluster_ids.size(); }
    inline uint32_t getTotalPoints() const { return n_total_points; }

    inline PointRef getPointRef(uint32_t idx)
    {
        return PointRef{&coords[idx * n_dim], &cluster_ids[idx], idx, n_dim};
    }
};