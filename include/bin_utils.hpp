#pragma once

#include "PointData.hpp"
#include "Centroids.hpp"

#include <fstream>

/// @brief Extracs all of the binary point data from the file pointed to by filename
/// @param filename 
/// @return 
PointData indvReadData(const std::string &filename);

/// @brief Extracts the binary point data pertaining to it's associated MPI Rank for all ranks.
/// @param filename 
/// @param n_mpi_procs 
/// @param mpi_rank 
/// @return 
PointData collectiveReadData(const std::string& filename, int n_mpi_procs, int mpi_rank);

void storeDataHeader(std::ofstream& file_stream, PointData data, uint32_t n_clusters);

void storeIterationData(std::ofstream& file_stream, PointData data, Centroids centroids);

void storeData(const std::string& filename, PointData data, Centroids centroids);

inline bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

/// @brief 
// Calculates the number of elements associated with a certain group, given the total number of elements between all groups.
/// @param n_total_points 
/// @param n_ranks 
/// @param rank 
/// @return 
uint32_t getGroupSize(uint32_t n_total_points, int n_ranks, int rank);
/// @brief 
// Calculates the offset (in element count) to the first element associated with a group, where all elements are divied evenly
// between groups. see ::getGroupSize
/// @param n_total_points 
/// @param n_ranks 
/// @param rank 
/// @return 
uint32_t getGroupOffset(uint32_t n_total_points, int n_ranks, int rank);


uint32_t getGroupForIndex(uint32_t index, uint32_t n_total_points, int n_ranks);