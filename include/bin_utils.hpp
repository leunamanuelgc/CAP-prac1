#pragma once

#include "Cluster.hpp"
#include "structs.hpp"

#include <fstream>

/**
 * Read data from the file: "salida" and stores it inside a dataResult struct.
 * @returns dataResult
 */
PointData readData(const std::string& filename, int n_mpi_procs, int mpi_rank);

void storeDataHeader(std::ofstream& file_stream, PointData data, std::vector<Cluster> clusters);

void storeIterationData(std::ofstream& file_stream, PointData data, std::vector<Cluster> clusters);

void storeData(const std::string& filename, PointData data, std::vector<Cluster> clusters);

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