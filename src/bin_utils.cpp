#include "bin_utils.hpp"
#include "Centroids.hpp"
#include "log_utils.hpp"

using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

//#define PRINT

/// @brief Single process reads entire input binary file.
/// @param filename 
/// @return 
PointData indvReadData(const string &filename)
{
    
    ifstream file(filename, ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    uint32_t n_rows, n_cols;
    // Leer el número de filas y columnas
    file.read(reinterpret_cast<char*>(&n_rows), sizeof(n_rows));
    file.read(reinterpret_cast<char*>(&n_cols), sizeof(n_cols));
    
    PointData data(n_rows,n_rows,n_cols);

    // Leer y guardar los puntos del fichero binario
    for (int i = 0; i < n_rows; i++) // each row represents a point's coordinates
    {
        file.read(reinterpret_cast<char*>(data.getPointCoords(i)), n_cols * sizeof(float));
    }

    return data;
};

PointData collectiveReadData(const string &filename, int n_mpi_procs, int mpi_rank)
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

    #ifdef PRINT
    std::cout << "File '" << filename << "' opened." << std::endl;
    #endif

    if (mpi_rank == 0)
    {
        MPI_File_read(fh, header, 2, MPI_UINT32_T, &frs);
    }
    MPI_Bcast(header, 2, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    #ifdef PRINT
    std::cout << "PointData header broadcasted successfully (points:" << header[0] << ", dim:" << header[1] << ")" << std::endl;
    #endif

    uint32_t n_rows=header[0], n_cols=header[1];
    auto n_local_p = getGroupSize(n_rows, n_mpi_procs, mpi_rank);
    auto local_p_offset = getGroupOffset(n_rows, n_mpi_procs, mpi_rank);

    PointData data(n_local_p, n_rows, n_cols);
    MPI_Offset point_coords_offset = sizeof(header) + local_p_offset * n_cols * sizeof(float);
    MPI_File_read_at(fh, point_coords_offset, data.coords.data(), n_local_p * n_cols, MPI_FLOAT, &frs);
    
    MPI_File_close(&fh);

    #ifdef PRINT
    std::cout << "Read point data (" << data.cluster_ids.size() << " points, " << data.coords.size() << " values) successfully!" << std::endl;
    #endif

    return data;
};

void storeDataHeader(ofstream &file_stream, PointData data, uint32_t n_clusters)
{
    file_stream.write(reinterpret_cast<char *>(&n_clusters), sizeof(n_clusters));
    auto np = data.getNPoints();
    file_stream.write(reinterpret_cast<char *>(&np), sizeof(np));
    auto dim = data.getDim();
    file_stream.write(reinterpret_cast<char *>(&dim), sizeof(dim));
}

void storeIterationData(ofstream &file_stream, PointData data, Centroids centroids)
{
    uint32_t n_clusters = centroids.n;
    // Store cluster centroid and id data
    for (int i = 0; i < n_clusters; i++)
    {
        file_stream.write(reinterpret_cast<char *>(&centroids[i][0]), sizeof(float) * data.getDim());
        auto cluster_id = i;
        file_stream.write(reinterpret_cast<char *>(&cluster_id), sizeof(cluster_id));
    }
    // Store point coordinates and cluster_id data
    for (int i = 0; i < data.getNPoints(); i++)
    {
        file_stream.write(reinterpret_cast<char *>(data.getPointCoords(i)), sizeof(float) * data.getDim());
        file_stream.write(reinterpret_cast<char *>(data.getPointClusterID(i)), sizeof(int));
    }
}

void storeData(const string &filename, PointData data, Centroids centroids)
{
    ofstream file(filename, ios::binary);
    if (!file)
    {
        throw std::runtime_error("Could not open file: " + filename);
    }
    storeDataHeader(file, data, centroids.n);
    storeIterationData(file, data, centroids);
}

uint32_t getGroupSize(uint32_t n_elements, int n_groups, int group)
{
    uint32_t base_amount = n_elements / n_groups;
    uint32_t remainder = n_elements % n_groups;
    return base_amount + (group < remainder ? 1 : 0);
}

uint32_t getGroupOffset(uint32_t n_elements, int n_groups, int group)
{
    int base_amount = n_elements / n_groups;
    int remainder = n_elements % n_groups;
    return group * base_amount + std::min(group, remainder);
}

uint32_t getGroupForIndex(uint32_t index, uint32_t n_elements, int n_groups)
{
    uint32_t base_amount = n_elements / n_groups;
    uint32_t remainder = n_elements % n_groups;

    // The first `remainder` groups have size (base_amount + 1)
    // Their offsets are 0, (base_amount+1), 2*(base_amount+1), ...
    if (index < (base_amount + 1) * remainder) {
        return index / (base_amount + 1);
    } else {
        return (index - ((base_amount + 1) * remainder)) / base_amount + remainder;
    }
}
