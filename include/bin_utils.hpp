#pragma once

#include "Cluster.hpp"
#include "structs.hpp"

#include <fstream>

/**
 * Read data from the file: "salida" and stores it inside a dataResult struct.
 * @returns dataResult
 */
PointData readData(const std::string& filename);

void storeDataHeader(std::ofstream& file_stream, PointData data, std::vector<Cluster> clusters);

void storeIterationData(std::ofstream& file_stream, PointData data, std::vector<Cluster> clusters);

void storeData(const std::string& filename, PointData data, std::vector<Cluster> clusters);

inline bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}