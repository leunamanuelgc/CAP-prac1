#include "bin_utils.hpp"

#include "log_utils.hpp"

using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

PointData readData(const string &filename, int n_mpi_procs, int mpi_rank)
{
    uint32_t header[2];
    MPI_File fh;
    MPI_Status frs;

    if (int err = MPI_File_open(MPI_COMM_WORLD, filename.data(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh); err != MPI_SUCCESS)
    {
        std::cout << "Could not open file: " << filename << std::endl;
        char err_str[MPI_MAX_ERROR_STRING];
        int len;
        MPI_Error_string(err, err_str, &len);
        std::cout << "MPI ERROR: " << std::string(err_str, len) << endl;
        
        MPI_Abort(MPI_COMM_WORLD, err);
    }
    if (mpi_rank == 0)
    {
        MPI_File_read(fh, header, 2, MPI_UINT32_T, &frs);
    }
    MPI_Bcast(header, 2, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    uint32_t n_rows=header[0], n_cols=header[1];
    auto n_local_p = getGroupSize(n_rows, n_mpi_procs, mpi_rank);
    auto local_p_offset = getGroupOffset(n_rows, n_mpi_procs, mpi_rank);
    
    vector<float> floats_v(n_local_p * n_cols);
    MPI_Offset point_coords_offset = sizeof(header) + local_p_offset * n_cols * sizeof(float);
    MPI_File_read_at(fh, point_coords_offset, floats_v.data(), n_local_p * n_cols, MPI_FLOAT, &frs);
    // DEBUG
    /*if (mpi_rank == -1)
    {
        vector<float> all_values(n_cols * n_rows);
        MPI_File_read_at(fh, point_coords_offset, all_values.data(), n_rows * n_cols, MPI_FLOAT, &frs);
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "All vectors(" << all_values.size()/n_cols << " floats) inside input file: " << std::endl;
        printVectors(all_values, n_cols);
    }
    else
    {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Local vectors(" << floats_v.size() << " floats) read from input file: " << std::endl;
        printVectors(floats_v, n_cols);
    }*/
    MPI_File_close(&fh);

    PointData data(n_local_p, n_rows, n_cols);
    //std::cout << mpi_rank << "\t N LOCAL Points: " << n_local_p << "\t First Point IDX: " << local_p_offset << std::endl;
    // Leer y guardar los puntos del fichero binario
    #pragma omp parallel for
    for (uint32_t i = 0; i < n_local_p; i++) // each row represents a point's coordinates
    {
        data.points[i].setID(i+local_p_offset);
        data.points[i].copyValueMemData(floats_v, i*n_cols, n_cols);
    }

    std::cout << "Read point data (" << data.points.size() << " points, " << floats_v.size() << " values) successfully!" << std::endl;

    return data;
};

void storeDataHeader(ofstream &file_stream, PointData data, vector<Cluster> clusters)
{
    uint32_t n_clusters = clusters.size();
    file_stream.write(reinterpret_cast<char *>(&n_clusters), sizeof(n_clusters));
    file_stream.write(reinterpret_cast<char *>(&data.n_total_points), sizeof(data.n_total_points));
    file_stream.write(reinterpret_cast<char *>(&data.n_dim), sizeof(data.n_dim));
}

void storeIterationData(ofstream &file_stream, PointData data, vector<Cluster> clusters)
{
    uint32_t n_clusters = clusters.size();
    // Store cluster centroid and id data
    for (int i = 0; i < n_clusters; i++)
    {
        file_stream.write(reinterpret_cast<char *>(&clusters[i].getCentroid()[0]), sizeof(float) * clusters[i].getDim());
        auto cluster_id = clusters[i].getID();
        file_stream.write(reinterpret_cast<char *>(&cluster_id), sizeof(clusters[i].getID()));
    }
    // Store point coordinates and cluster_id data
    for (int i = 0; i < data.n_total_points; i++)
    {
        file_stream.write(reinterpret_cast<char *>(&data.points[i].getValues()[0]), sizeof(float) * data.n_dim);
        auto cluster_id = data.points[i].getClusterID();
        file_stream.write(reinterpret_cast<char *>(&cluster_id), sizeof(cluster_id));
    }
}

void storeData(const string &filename, PointData data, vector<Cluster> clusters)
{
    ofstream file(filename, ios::binary);
    if (!file)
    {
        throw std::runtime_error("Could not open file: " + filename);
    }
    storeDataHeader(file, data, clusters);
    storeIterationData(file, data, clusters);
}

uint32_t getGroupSize(uint32_t n_total_points, int n_ranks, int rank)
{
    uint32_t base_amount = n_total_points / n_ranks;
    uint32_t remainder = n_total_points % n_ranks;
    return base_amount + (rank < remainder ? 1 : 0);
}

uint32_t getGroupOffset(uint32_t n_total_points, int n_ranks, int rank)
{
    int base_amount = n_total_points / n_ranks;
    int remainder = n_total_points % n_ranks;
    return rank * base_amount + std::min(rank, remainder);
}
