#pragma once

#include "Point.hpp"

#include <cstdint>
#include <vector>
#include <mpi.h>
#include <memory>

// Place smaller fields towards the end
// just in case. For lower padding and
// better memory performance

struct PointData
{
    std::vector<Point> points;
    uint32_t n_total_points;
    uint32_t n_dim;

    inline PointData() {}
    inline PointData(uint32_t n_points, uint32_t n_total_points, uint32_t n_dim) :
        points(n_points, Point(n_dim)),
        n_total_points(n_total_points),
        n_dim(n_dim) {}
};

struct CentroidDiff
{
    std::vector<float> add_points_sum;
    std::vector<float> rem_points_sum;
    uint32_t add_points_count;
    uint32_t rem_points_count;

    // inline CentroidDiff() {}
    inline CentroidDiff(uint32_t dim) : add_points_sum(dim, 0),
                                        rem_points_sum(dim, 0),
                                        add_points_count(0),
                                        rem_points_count(0) {}

    inline size_t nFlatBytes() const { return add_points_sum.size() * sizeof(float) * 2 + sizeof(uint32_t) * 2; }
    inline static uint32_t dim(size_t n_flat_bytes) { return (n_flat_bytes - 2 * sizeof(uint32_t)) / (2 * sizeof(float)); }
    // Function to flatten data into a byte buffer
    void copyBytesIntoFlatBuff(std::byte *data_ptr) const;
    void copyFromFlatBytes(const std::byte *flat_byte_ptr);
    static void addFlatDiffs(void *a, void *b, int dim);
};

// Sum operator used in the reduction of the centroid diffs
void centroid_diff_sum_function(void *inputBuffer, void *outputBuffer, int *len, MPI_Datatype *datatype);
