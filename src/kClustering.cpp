#include "PointData.hpp"
#include "PointRef.hpp"
#include "Centroids.hpp"
#include "CentroidDiffs.hpp"
#include "CentroidDiffRef.hpp"
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
void assignPointToCluster(PointRef &p, int new_cluster_id, CentroidDiffs &centroid_diffs);
int shareAndApplyCentroidDiffs(int k, Centroids &centroids, CentroidDiffs &centroid_diffs, std::vector<uint32_t> &cluster_point_counts);

//#define DEBUG
#define VERBOSE 2
#define LOG 1

//#define PRINT LOG
#define PRINT VERBOSE

#define INPUT_DATA "./build/data/salida"
#define OUTPUT_DATA "./build/data/cluster_data"

// Bug: high point counts release builds might segfault if store iterations is activated
//      but storeiterations shouldn't be used for high point counts anyways, since it's
//      main purpose is to use in datavis to check if the algorithm is running correctly.
//#define STORE_ITERATIONS

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    const int N_CLUSTERS = 5;

    int n_mpi_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    #ifdef DEBUG
    std::cout << "Rank: " << mpi_rank << "/" << n_mpi_procs-1 << std::endl;
    #endif

    // Obtencion de los puntos por cada proceso
    PointData data = collectiveReadData(INPUT_DATA, n_mpi_procs, mpi_rank);
    
    // Hacer backup de los datos de salida preexistentes
    if (mpi_rank == 0 && fileExists(OUTPUT_DATA) && argv[0] == "b")
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

void kMeans(PointData data, int k, int n_mpi_procs, int mpi_rank)
{
    #if PRINT >= VERBOSE
    std::cout << " Running distributed K-Means algorithm..." << std::endl;
    #endif
    
    const int MAX_ITERATIONS = 100;
    const float MIN_PTS_MVMT_PCT = 0.05;

    // Benchmarking
    auto start = chrono::high_resolution_clock::now();
    auto end = start;
    auto total_start = start;
    auto total_end = end;
    chrono::duration<double> elapsed, total_elapsed;
    double elapsed_time;
    double iter_t;
    double total_t;

    // Representa el n de puntos correspondiente al proceso mpi
    const int n_local_p = data.getNPoints();
    const int p_offset = getGroupOffset(data.getTotalPoints(), n_mpi_procs, mpi_rank);
    const auto dim = data.getDim();
    const auto n_total_points = data.getTotalPoints();

    
#ifdef STORE_ITERATIONS
    const int STORE_RANK = 0;
    PointData merged_data;
    ofstream output_file;
    if (mpi_rank == STORE_RANK)
    {
        merged_data = indvReadData(INPUT_DATA);
        output_file = ofstream(OUTPUT_DATA, ios::binary);
        if (!output_file) {
            throw std::runtime_error("Could not open file: " + (string)OUTPUT_DATA);
        }
        // And write header
        storeDataHeader(output_file, merged_data, k);
    }
#endif

    int iteration = 0;
    bool convergence_criterion = false;

    // Centroid representing the clusters' "average point" or "middle".
    // Ought to be shared and correctly synchronized between processes
    Centroids centroids(k, dim);

    vector<uint32_t> cluster_point_counts(k, 0);

    // Local accumulators for each processes' changes in cluster centroids (both local and foreign)
    CentroidDiffs centroid_diffs(k, dim);

    // Inicializar los clusters (centroides)
    for (int i = 0; i < n_local_p; i++)
    {
        int k_init_idx = getGroupForIndex(i+p_offset, n_total_points, k);
        auto p = data.getPointRef(i);
        assignPointToCluster(p, k_init_idx, centroid_diffs);
    }
    
    shareAndApplyCentroidDiffs(k, centroids, centroid_diffs, cluster_point_counts);

    
    // DEBUG
    #ifdef DEBUG
    //if (mpi_rank == 0)
    {
        for (int i = 0; i < k; i++)
        {
            std::cout << "CENTROID_DIFF[" << i << "]: -(" << centroid_diffs[i].getRemPointsCount() << "):";
            printVector(centroid_diffs[i].rem_points_sum, dim);
            std::cout << "CENTROID_DIFF[" << i << "]: +(" << centroid_diffs[i].getAddPointsCount() << "):";
            printVector(centroid_diffs[i].add_points_sum, dim);
            std::cout << endl;
        }

        for (int i = 0; i < k; i++)
        {
            std::cout << "CENTROID[" << i << "]: \t";
            printVector(centroids[i], dim);
            std::cout << std::endl;
        }
    }
    #endif

    while (iteration < MAX_ITERATIONS && !convergence_criterion)
    {
        auto it_start = chrono::high_resolution_clock::now();
        auto it_end = it_start;

        double dist_calc_t = 0;
        centroid_diffs = CentroidDiffs(k, dim);

        // Should probably parallelize THIS outer for loop with OpenMP rather than the inner levels,
        // as it usually provides better results the further out (to an extent) the parallelized loop is.
        // NOTE: though it would call for  a couple reductions, it's probably worth the hassle.
        for (int i = 0; i < n_local_p; i++)
        {
            PointRef p = data.getPointRef(i);
            double min_dist = numeric_limits<double>::infinity();
            int closest_cluster_id;

            start = chrono::high_resolution_clock::now();

            // Calcular las distancias entre los puntos y los centroides
            // *emabarrassingly parallel
            #pragma omp parallel shared(p, centroids, min_dist)
            {
                float th_min_dist = numeric_limits<float>::infinity();
                int th_closest_cluster_id;

                #pragma omp for nowait
                for (int j = 0; j < k; j++)
                {
                    float dist = sqrDist(p.getValues(), centroids[j], dim);
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

            // Assign point to the new closest cluster
            if (p.getClusterID() != closest_cluster_id)
            {
                assignPointToCluster(p, closest_cluster_id, centroid_diffs);
            }
        }

        int moved_points = shareAndApplyCentroidDiffs(k, centroids, centroid_diffs, cluster_point_counts);

        iteration++;

        if (moved_points < (n_total_points * MIN_PTS_MVMT_PCT))
        {
            convergence_criterion = true;
        }

        it_end = chrono::high_resolution_clock::now();

        elapsed = it_end - it_start;
        elapsed_time = elapsed.count();

        iter_t = elapsed_time;

    // PRINTING
    #ifdef DEBUG
        if (mpi_rank == 0)
        {
            printCentroidsData(centroids, centroid_diffs, cluster_point_counts);
        }
    #endif
    #if PRINT >= VERBOSE
        cout << "ITERATION " << iteration - 1 << endl;
        printLocalPointInfo(data, iteration-1,  n_mpi_procs, mpi_rank);
        printMovedPoints(moved_points, MIN_PTS_MVMT_PCT, iteration, convergence_criterion);
        cout << endl;
    #endif
    #if PRINT >= LOG
        printBenchmarkCSV(iteration - 1, iter_t, dist_calc_t);
    #endif

    #ifdef STORE_ITERATIONS
        if (mpi_rank != STORE_RANK)
        {
            MPI_Gatherv(data.cluster_ids.data(), data.getNPoints(), MPI_INT, NULL, NULL, NULL, MPI_INT, STORE_RANK, MPI_COMM_WORLD);
        } else {
            int recvcounts[k];
            int displacements[k];
            #pragma omp parallel for
            for (int i = 0; i < k; i++)
            {
                recvcounts[i] = getGroupSize(n_total_points, n_mpi_procs, i);
                displacements[i] = getGroupOffset(n_total_points, n_mpi_procs, i);
            }
            MPI_Gatherv(data.cluster_ids.data(), data.getNPoints(), MPI_INT, merged_data.cluster_ids.data(), recvcounts, displacements, MPI_INT, STORE_RANK, MPI_COMM_WORLD);
            storeIterationData(output_file, merged_data, centroids);
        }
    #endif
    }

#if PRINT >= LOG
    total_end = chrono::high_resolution_clock::now();
    total_elapsed = total_end - total_start;
    total_t = total_elapsed.count();
    std::cout << "total run t., " << total_t
    << ", n. iters., " << iteration
    << ", n. points, " << data.getTotalPoints()
    << ", n. clusters, " << k
    << std::endl;
#endif


}

void assignPointToCluster(PointRef &p, int new_cluster_id, CentroidDiffs &centroid_diffs)
{
    //std::cout << "assigning point(" << p.getID() << ") to cluster(" << new_cluster_id << ")\t old cluster: " << p.getClusterID() << std::endl;
    int old_cluster_id = p.getClusterID();
    // Remove from previous cluster
    if (old_cluster_id != NONE_CLUSTER)
    {
        // Centroid diff calculation
        centroid_diffs[old_cluster_id].incrementRemCount();
        addInto(centroid_diffs[old_cluster_id].rem_points_sum, p.getValues(), p.getDim());

    }
    // Cendtroid diff calculation
    centroid_diffs[new_cluster_id].incrementAddCount();
    addInto(centroid_diffs[new_cluster_id].add_points_sum, p.getValues(), p.getDim());

    p.setClusterID(new_cluster_id);
}

int shareAndApplyCentroidDiffs(int k, Centroids &centroids, CentroidDiffs &centroid_diffs, vector<uint32_t> &cluster_point_counts)
{
    int moved_points = 0;

    // Fully exchange centroid diff and point movement information
    MPI_Allreduce(MPI_IN_PLACE, centroid_diffs.add_rem_points_sums.data(), centroid_diffs.add_rem_points_sums.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, centroid_diffs.add_rem_points_counts.data(), centroid_diffs.add_rem_points_counts.size(), MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);

    // Readjust new cluster centroid (cluster's points' average)
    // #pragma omp parallel for default(private) <-- THIS WOULD CAUSE LINKING ERRORS:
    //                                              because of unexplicit variable initialization behaviour,
    //                                              as a result templated classes like those of std (i.e vector)
    //                                              will fail to produce default-constructors during compilation
    //                                              and at the linking step won't be able to find the corresponding
    //                                              constructors demanded by each thread's implicit default initialization
    //                                              of the privated variables.
    #pragma omp parallel for default(none) \
        shared(centroid_diffs, k, centroids, cluster_point_counts) \
        reduction(+ : moved_points)
    for (int i = 0; i < k; i++)
    {
        uint32_t prev_point_count = cluster_point_counts[i];
        cluster_point_counts[i] += centroid_diffs[i].getAddPointsCount() - centroid_diffs[i].getRemPointsCount();
        uint32_t new_point_count = cluster_point_counts[i];

        moved_points += centroid_diffs[i].getAddPointsCount();

        if (new_point_count > 0)
        {
            // PROBAR CAMBIAR PRAGMA AQUI A VER SI VA MEJOR
            for (int j = 0; j < centroids.dim; j++)
            {
                float add_points_mean_val = 0;
                float rem_points_mean_val = 0;
    
    
                if (centroid_diffs[i].getAddPointsCount() > (uint32_t)0)
                {
                    add_points_mean_val = centroid_diffs[i].add_points_sum[j] / centroid_diffs[i].getAddPointsCount();
                }
    
                if (centroid_diffs[i].getRemPointsCount() > (uint32_t)0)
                {
                    rem_points_mean_val = centroid_diffs[i].rem_points_sum[j] / centroid_diffs[i].getRemPointsCount();
                }
    
                centroids[i][j] =
                (
                    centroids[i][j] * prev_point_count
                    + add_points_mean_val * (centroid_diffs[i].getAddPointsCount())
                    - rem_points_mean_val * (centroid_diffs[i].getRemPointsCount())
                ) / new_point_count;
            }
        }
    }

    return moved_points;
}
