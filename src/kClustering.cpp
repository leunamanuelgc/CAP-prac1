#include "Point.hpp"
#include "Cluster.hpp"
#include "bin_utils.hpp"
#include "log_utils.hpp"

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>
#include <iomanip>
#include <chrono>
#include <array>
#include <omp.h>
#include <mpi.h>
#include <string.h>

void kMeans(PointData data, int k, int n_mpi_procs, int mpi_rank);

#define VERBOSE 2
#define LOG 1

// #define PRINT LOG
#define PRINT VERBOSE

#define INPUT_DATA "../data/salida"
#define OUTPUT_DATA "../data/cluster_data"

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    const int N_CLUSTERS = 10;

    int n_mpi_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    std::cout << "Rank: " << mpi_rank << "/" << n_mpi_procs-1 << std::endl;

    // Obtencion de los puntos por cada proceso
    PointData data = readData(INPUT_DATA, n_mpi_procs, mpi_rank);

    // Hacer backup de los datos de salida preexistentes
    if (fileExists(OUTPUT_DATA) && argv[0] == "b")
    {
        int backup_id = 0;
        bool file_created = false;
        while (!file_created)
        {
            std::string backup_name = OUTPUT_DATA + (std::string) "_" + std::to_string(backup_id);
            if (!fileExists(backup_name))
            {
                file_created = (rename(OUTPUT_DATA, backup_name.c_str()) == 0);
            }
            backup_id++;
        }
    }

    srand(time(NULL));

    kMeans(data, N_CLUSTERS, n_mpi_procs, mpi_rank);

    MPI_Finalize();

    return 0;
}

using namespace std;

int nLocalK(int k, int n_mpi_procs, int rank)
{
    const int last_mpi_rank = n_mpi_procs - 1;
    return (int)k / n_mpi_procs + rank == last_mpi_rank ? k % n_mpi_procs : 0;
}

int firstLocalK(int k, int n_mpi_procs, int rank)
{
    return k / n_mpi_procs * rank;
}

void kMeans(PointData data, int k, int n_mpi_procs, int mpi_rank)
{
    const int MAX_ITERATIONS = 2000;
    const float MIN_PTS_MVMT_PCT = 0.05;
    // Benchmarking
    auto start = chrono::high_resolution_clock::now();
    auto end = start;
    chrono::duration<double> elapsed;
    double elapsed_time;
    double total_dist_calc_t;
    double total_centroid_calc_t;
    double iter_t;

    // Representa el n de local_clusters correspondiente al proceso mpi
    //  * give remainder to the last to simplify indexing to the first for each process
    const int n_local_k = nLocalK(k, n_mpi_procs, mpi_rank);
    const int first_local_k = firstLocalK(k, n_mpi_procs, mpi_rank);

    // Representa el n de puntos correspondiente al proceso mpi
    const int n_local_p = data.points.size();

    /// TODO: Binary output data file writing and formatting
    /*
    // Open file for iteration output on MPI rank 0
    if (mpi_rank == 0)
    {
        ofstream file("clustered_data", ios::binary);
        if (!file) {
            throw std::runtime_error("Could not open file: clustered_data");
        }
        // And write header
        storeDataHeader(file, data, local_clusters);
    }
    */

    int iteration = 0;
    bool convergence_criterion = false;

    // Centroid representing the clusters' "average point" or "middle".
    // Ought to be shared and correctly synchronized between processes
    vector<vector<float>> centroids(k, vector<float>(data.n_dim));
    vector<uint32_t> cluster_point_counts(k, 0);

    // Local accumulators for each processes' changes in cluster centroids (both local and foreign)
    vector<CentroidDiff> centroidDiffs(k, CentroidDiff(data.n_dim));

    // Register reduction operation for centroid diffs to mpi
    // (sums the corresponding CentroidDiffs individually)
    MPI_Op op_sum_centroid_diffs;
    MPI_Op_create(centroid_diff_sum_function, true, &op_sum_centroid_diffs);

    // Register the flat centroid diff datatype to mpi
    MPI_Datatype t_flat_bytes_centroid_diff;
    size_t n_flat_bytes = CentroidDiff::nFlatBytes(data.n_dim);
    MPI_Type_contiguous(n_flat_bytes, MPI_BYTE, &t_flat_bytes_centroid_diff);
    MPI_Type_commit(&t_flat_bytes_centroid_diff);
    // working centroid diffs flattened array
    vector<byte> flat_centroid_diffs(CentroidDiff::nFlatBytes(data.n_dim) * k);

    // Inicializar los clusters (centroides)
    {
        // Contains the coordinates of each centroid, ordered by cluster id.
        vector<float> bcst_centroids(data.n_dim * k);
        vector<int> n_local_ks(k);
        vector<int> first_local_ks(k, 0);

#pragma omp parallel for
        for (int i = 0; i < k; i++)
        {
            n_local_ks[i] = nLocalK(k, n_mpi_procs, i);
            first_local_ks[i] = firstLocalK(k, n_mpi_procs, i);

            if (i >= first_local_k && i < (first_local_k + n_local_k))
            {
                // Pick out equally spaced points from the data input vector
                // (i.e random points, except not bc that's actually how input data is structured in the clustered data generator)
                int centroid_idx = (n_local_p / n_local_k) * i + n_local_p / n_local_k / 2;

                // RANDOM INIT
                // int centroid_idx = i;

                // RANDOM INIT
                // int centroid_idx = (int)((float)rand()/RAND_MAX * n_local_p) + first_local_p;
                // cout << "Centroid[" << i << "]: \t" << centroid_idx << endl;
                auto centroid = data.points[centroid_idx].getValues();
                // lo suyo probablemente seria usar std::copy, pero
                // no me apetece comerme complejidad O(N) por la cara.
                memcpy(&bcst_centroids[i * data.n_dim], centroid.data(), data.n_dim * sizeof(float));
            }
        }

        // Esto lo deberiamos quitar y hacer que lo pueda calcular de forma determinista y ya pero bueno...
        // Communicate and sync centroid initialization values
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, bcst_centroids.data(), n_local_ks.data(), first_local_ks.data(), MPI_FLOAT, MPI_COMM_WORLD);

#pragma omp parallel for
        for (int i = 0; i < k; i++)
        {
            centroids[i].assign(bcst_centroids[i * data.n_dim], bcst_centroids[(i + 1) * data.n_dim]);
        }
    }

    while (iteration < MAX_ITERATIONS && !convergence_criterion)
    {
        auto it_start = chrono::high_resolution_clock::now();
        auto it_end = it_start;
        int moved_points = 0;

        double dist_calc_t = 0;
        centroidDiffs = vector<CentroidDiff>(k, CentroidDiff(data.n_dim));

        // Should probably parallelize THIS outer for loop with OpenMP rather than the inner levels,
        // as it usually provides better results the further out (to an extent) the parallelized loop is.
        // NOTE: though it would call for  a couple reductions, it's probably worth the hassle.
        for (int i = 0; i < n_local_p; i++)
        {
            Point &p = data.points[i];
            double min_dist = numeric_limits<double>::infinity();
            int closest_cluster_id;

            start = chrono::high_resolution_clock::now();

// Calcular las distancias entre los puntos y los centroides
// *emabarrassingly parallel
#pragma omp parallel shared(p, centroids, min_dist)
            {
                double th_min_dist = numeric_limits<double>::infinity();
                int th_closest_cluster_id;

#pragma omp for nowait
                for (int j = 0; j < k; j++)
                {
                    double dist = sqrDist(p.getValues(), centroids[j]);
                    if (dist < th_min_dist)
                    {
                        th_min_dist = dist;
                        th_closest_cluster_id = j;
                    }
                }

// Bernat:
// critical section extracted from thread group loop
// to avoid reocurring locks between threads, since
// threads are handed out a portion of the loop to
// run, the critical section would be hit constantly
// greatly increasing wait times for the 'mutex' to clear
#pragma omp critical
                {
                    // Comprobar cual es la menor distancia entre un punto y los centroides
                    if (th_min_dist < min_dist)
                    {
                        min_dist = th_min_dist;
                        closest_cluster_id = th_closest_cluster_id;
                    }
                }
            }
            end = chrono::high_resolution_clock::now();
            elapsed = end - start;
            elapsed_time = elapsed.count();
            dist_calc_t += elapsed_time;

            //  sync

            // Assign point to the new closest cluster
            if (p.getClusterID() != closest_cluster_id)
            {
                int old_cluster_id = p.getClusterID();
                // Remove from previous cluster
                if (old_cluster_id != NONE_CLUSTER)
                {
                    // clusters[old_cluster_id].removePointByID(p.getID());

                    // Centroid diff calculation
                    centroidDiffs[old_cluster_id].rem_points_count++;
                    for (int j = 0; j < centroidDiffs[old_cluster_id].rem_points_sum.size(); j++)
                    {
                        centroidDiffs[old_cluster_id].rem_points_sum[j] += p.getValueAt(j);
                    }
                }

                // clusters[closest_cluster_id].addPoint(p);

                // Cendtroid diff calculation
                centroidDiffs[closest_cluster_id].add_points_count++;
                for (int j = 0; j < centroidDiffs[closest_cluster_id].add_points_sum.size(); j++)
                {
                    centroidDiffs[closest_cluster_id].add_points_sum[j] += p.getValueAt(j);
                }

                p.setClusterID(closest_cluster_id);
            }
        }

#pragma omp parallel for default(none) shared(k, flat_centroid_diffs, centroidDiffs, n_flat_bytes)
        for (int i = 0; i < k; i++)
        {
            byte *fcd_ptr = &flat_centroid_diffs[i * n_flat_bytes];
            centroidDiffs[i].copyBytesIntoFlatBuff(fcd_ptr);
        }

        // Fully exchange centroid diff and point movement information
        MPI_Allreduce(MPI_IN_PLACE, flat_centroid_diffs.data(), k, t_flat_bytes_centroid_diff, op_sum_centroid_diffs, MPI_COMM_WORLD);

// Readjust new cluster centroid (cluster's points' average)
// #pragma omp parallel for default(private) <-- THIS WOULD CAUSE LINKING ERRORS:
//                                              because of unexplicit variable initialization behaviour,
//                                              as a result templated classes like those of std (i.e vector)
//                                              will fail to produce default-constructors during compilation
//                                              and at the linking step won't be able to find the corresponding
//                                              constructors demanded by each thread's implicit default initialization
//                                              of the privated variables.
#pragma omp parallel for default(none)                                                           \
    shared(centroidDiffs, flat_centroid_diffs, n_flat_bytes, k, centroids, cluster_point_counts) \
    reduction(+ : moved_points)
        for (int i = 0; i < k; i++)
        {
            centroidDiffs[i].copyFromFlatBytes(&flat_centroid_diffs[i * n_flat_bytes]);

            uint32_t prev_point_count = cluster_point_counts[i];
            cluster_point_counts[i] += centroidDiffs[i].add_points_count - centroidDiffs[i].rem_points_count;
            uint32_t new_point_count = cluster_point_counts[i];

            moved_points += centroidDiffs[i].add_points_count;

            if (centroidDiffs[i].add_points_count + centroidDiffs[i].rem_points_count == 0)
                continue;

            // copy old centroid
            vector<float> new_centroid(centroids[i]);
            for (int j = 0; j < new_centroid.size(); j++)
            {
                centroids[i][j] *= (float)prev_point_count / new_point_count;
            }

            if (centroidDiffs[i].add_points_count > 0)
            {
                for (int j = 0; j < new_centroid.size(); j++)
                {
                    float add_points_mean_val = centroidDiffs[i].add_points_sum[j] / centroidDiffs[i].add_points_count;
                    centroids[i][j] +=
                        add_points_mean_val *
                        (float)centroidDiffs[i].add_points_count / new_point_count;
                }
            }
            if (centroidDiffs[i].rem_points_count > 0)
            {
                for (int j = 0; j < centroids[i].size(); j++)
                {
                    float rem_points_mean_val = centroidDiffs[i].rem_points_sum[j] / centroidDiffs[i].rem_points_count;
                    centroids[i][j] -=
                        rem_points_mean_val *
                        (float)centroidDiffs[i].rem_points_count / prev_point_count;
                }
            }
        }

        iteration++;

        if (moved_points < (data.n_total_points * MIN_PTS_MVMT_PCT))
        {
            convergence_criterion = true;
        }

        it_end = chrono::high_resolution_clock::now();

        elapsed = it_end - it_start;
        elapsed_time = elapsed.count();

        iter_t = elapsed_time;

// PRINTING
#if PRINT >= VERBOSE
        cout << "ITERATION " << iteration - 1 << endl;
        printLocalPointInfo(data.points, iteration-1,  n_mpi_procs, mpi_rank);
        printMovedPoints(moved_points, MIN_PTS_MVMT_PCT, iteration, convergence_criterion);
        cout << endl;
#elif PRINT >= LOG
        printBenchmarkCSV(iteration - 1, iter_t, dist_calc_t);
#endif

        // storeIterationData(file, data, local_clusters);
    }

// Trying out per-iteration data storage
// storeData("clustered_data", data, local_clusters);

// Benchmarking
#if PRINT >= VERBOSE
    cout << "Total distance compute time: " << total_dist_calc_t << "(avg. per iteration: " << total_dist_calc_t / (iteration - 1) << ")\n";
#elif PRINT >= LOG
// cout << "total dist. compute T" << ", " << "avg. T per iter" << endl;
// cout << total_dist_calc_t << ", " << total_dist_calc_t/(iteration-1) << endl;
#endif
}
