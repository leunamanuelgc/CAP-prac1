#include "bin_utils.hpp"

using namespace std;

PointData readData(const string& filename, int n_mpi_procs, int mpi_rank)
{
    ifstream file(filename, ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    
    PointData data;
    uint32_t n_rows, n_cols;
    
    // Leer el número de filas y columnas
    file.read(reinterpret_cast<char*>(&n_rows), sizeof(n_rows));
    file.read(reinterpret_cast<char*>(&n_cols), sizeof(n_cols));
    data.n_total_points = n_rows;
    data.n_dim = n_cols;
    
    const int last_mpi_rank = n_mpi_procs-1;
    // Representa el n de puntos correspondiente al proceso mpi
    const int n_local_p = (int)data.n_total_points/n_mpi_procs + mpi_rank == last_mpi_rank ? data.n_total_points % n_mpi_procs : 0;
    const int first_local_p = (int)data.n_total_points/n_mpi_procs * mpi_rank;

    const size_t read_bytes_offset = first_local_p * n_cols * sizeof(float);

    // Offset file reading to the start of the processes' first local point index,
    // starting from AFTER the "header" containing 'n_total_points' and 'n_dim'
    file.seekg(read_bytes_offset, ios::cur);

    // Leer y guardar los puntos del fichero binario
    for (int i = 0; i < n_local_p; i++) // each row represents a point's coordinates
    {
        vector<float> values(n_cols);
        file.read(reinterpret_cast<char*>(&values[0]), n_cols * sizeof(float));
        Point p(i, NONE_CLUSTER, n_cols, values);
        data.points.push_back(move(p));
    }

    return data;
};

void storeDataHeader(ofstream& file_stream, PointData data, vector<Cluster> clusters)
{
    uint32_t n_clusters = clusters.size();
    file_stream.write(reinterpret_cast<char*>(&n_clusters), sizeof(n_clusters));
    file_stream.write(reinterpret_cast<char*>(&data.n_total_points), sizeof(data.n_total_points));
    file_stream.write(reinterpret_cast<char*>(&data.n_dim), sizeof(data.n_dim));
}

void storeIterationData(ofstream& file_stream, PointData data, vector<Cluster> clusters)
{
    uint32_t n_clusters = clusters.size();
    // Store cluster centroid and id data
    for (int i = 0; i < n_clusters; i++)
    {
        file_stream.write(reinterpret_cast<char*>(&clusters[i].getCentroid()[0]), sizeof(float) * clusters[i].getDim());
        auto cluster_id = clusters[i].getID();
        file_stream.write(reinterpret_cast<char*>(&cluster_id), sizeof(clusters[i].getID()));
    }
    // Store point coordinates and cluster_id data
    for (int i = 0; i < data.n_total_points; i++)
    {
        file_stream.write(reinterpret_cast<char*>(&data.points[i].getValues()[0]), sizeof(float) * data.n_dim);
        auto cluster_id = data.points[i].getClusterID();
        file_stream.write(reinterpret_cast<char*>(&cluster_id), sizeof(cluster_id));
    }
}

void storeData(const string& filename, PointData data, vector<Cluster> clusters)
{
    ofstream file(filename, ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    storeDataHeader(file, data, clusters);
    storeIterationData(file, data, clusters);
}
