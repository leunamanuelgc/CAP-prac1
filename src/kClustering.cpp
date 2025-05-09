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


void kMeans(PointData data, int k);

#define VERBOSE 2
#define LOG 1

#define PRINT LOG
//#define PRINT VERBOSE

#define INPUT_DATA "../data/salida"
#define OUTPUT_DATA "../data/cluster_data"

int main(int argc, char **argv)
{
    const int N_CLUSTERS = 10;

    // Obtencion de los puntos
    PointData data = readData(INPUT_DATA);

    // Hacer backup de los datos de salida preexistentes
    if (fileExists(OUTPUT_DATA) && argv[0] == "b")
    {
        int backup_id = 0;
        bool file_created = false;
        while (!file_created)
        {
            std::string backup_name = OUTPUT_DATA + (std::string)"_" + std::to_string(backup_id);
            if (!fileExists(backup_name))
            {
                file_created = (rename(OUTPUT_DATA, backup_name.c_str()) == 0);
            }
            backup_id++;
        }
    }

    srand(time(NULL));

    kMeans(data, N_CLUSTERS);

    return 0;
}

using namespace std;


void kMeans(PointData data, int k)
{
    const int MAX_ITERATIONS = 30;
    const float MIN_PTS_MVMT_PCT = 0.05;
    // Benchmarking
	auto start = chrono::high_resolution_clock::now();
	auto end = start;
	chrono::duration<double> elapsed;
	double elapsed_time;
    double total_dist_calc_t;
    double total_centroid_calc_t;
    double iter_t;
    
    std::vector<Cluster> clusters(k, Cluster(data.n_dim));

    // Inicializar los clusters
    for (int i = 0; i < k; i++)
    {
        // Pick out equally spaced points from the data input vector (i.e random points, except not bc yes)
        int centroid_idx = (data.n_points / k) * i + data.n_points/k/2;

        // RANDOM INIT
        //int centroid_idx = i;

        // RANDOM INIT
        //int centroid_idx = (int)((float)rand()/RAND_MAX * data.n_points);
        //cout << "Centroid[" << i << "]: \t" << centroid_idx << endl;
        clusters[i] = Cluster(i, data.points[centroid_idx].getValues());
    }

    // Open file for iteration output
    ofstream file("clustered_data", ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: clustered_data");
    }
    // And write header
    storeDataHeader(file, data, clusters);

    int iteration = 0;
    bool convergence_criterion = false;

    vector<uint32_t> old_cluster_p_count(k, 0);
    while (iteration < MAX_ITERATIONS && !convergence_criterion)
    {
        auto it_start = chrono::high_resolution_clock:: now();
        auto it_end = it_start;
        int moved_points = 0;

        double dist_calc_t  = 0;
        vector<CentroidDiff> centroidDiffs(k, CentroidDiff(data.n_dim));

        for (int i = 0; i < data.n_points; i++)
        {
            Point& p = data.points[i];
            double min_dist = numeric_limits<double>::infinity();
            int closest_cluster_id;

            start = chrono::high_resolution_clock::now();
            
            // Calcular las distancias entre los puntos y los centroides de los k nodos
            // *emabarrassingly parallel
            #pragma omp parallel shared(p,clusters, min_dist)
            {
                double th_min_dist = numeric_limits<double>::infinity();
                int th_closest_cluster_id;
            
                #pragma omp for nowait
                for (int j = 0; j < k; j++)
                {
                    double dist = sqrDist(p.getValues(), clusters[j].getCentroid());
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
                    clusters[old_cluster_id].removePointByID(p.getID());
                    
                    // Centroid diff calculation
                    centroidDiffs[old_cluster_id].rem_points_count++;
                    for (int j = 0; j < centroidDiffs[old_cluster_id].rem_points_sum.size(); j++)
                    {
                        centroidDiffs[old_cluster_id].rem_points_sum[j] += p.getValueAt(j);
                    }
                }
                
                clusters[closest_cluster_id].addPoint(p);

                // Cendtroid diff calculation
                centroidDiffs[closest_cluster_id].add_points_count++;
                for (int j = 0; j < centroidDiffs[closest_cluster_id].add_points_sum.size(); j++)
                {
                    centroidDiffs[closest_cluster_id].add_points_sum[j] += p.getValueAt(j);
                }

                p.setClusterID(closest_cluster_id);
                moved_points++;
            }
        }

        // sync -> reduction (~average) (actually we will do a weithed
        //         mean calculation as a reduction method for the centroid)

        // Readjust new cluster centroid (cluster's points' average)
        //#pragma omp parallel for default(private) <-- THIS WOULD CAUSE LINKING ERRORS:
        //                                              because of unexplicit variable initialization behaviour,
        //                                              as a result templated classes like those of std (i.e vector)
        //                                              will fail to produce default-constructors during compilation
        //                                              and at the linking step won't be able to find the corresponding
        //                                              constructors demanded by each thread's implicit default initialization
        //                                              of the privated variables.
        #pragma omp parallel for
        for (int i = 0; i < k; i++)
        {
            if (centroidDiffs[i].add_points_count + centroidDiffs[i].rem_points_count == 0)
                continue;

            // copy old centroid
            vector<float> new_centroid(clusters[i].getCentroid());
            // maybe faster if initizlization is done inside the loop too
            for (int j = 0; j < new_centroid.size(); j++)
            {
                new_centroid[j] *= (float)old_cluster_p_count[i] / clusters[i].getNumPoints();
            }
            
            if (centroidDiffs[i].add_points_count > 0)
            {
                for (int j = 0; j < new_centroid.size(); j++)
                {
                    float add_points_mean_val = centroidDiffs[i].add_points_sum[j] / centroidDiffs[i].add_points_count;
                    new_centroid[j] +=
                        add_points_mean_val *
                        (float)centroidDiffs[i].add_points_count/clusters[i].getNumPoints();
                }
            }
            if (centroidDiffs[i].rem_points_count > 0)
            {
                for (int j = 0; j < new_centroid.size(); j++)
                {
                    float rem_points_mean_val = centroidDiffs[i].rem_points_sum[j] / centroidDiffs[i].rem_points_count;
                    new_centroid[j] -=
                        rem_points_mean_val *
                        (float)centroidDiffs[i].rem_points_count/clusters[i].getNumPoints();
                }
            }
            
            clusters[i].setCentroid(new_centroid);
            old_cluster_p_count[i] = clusters[i].getNumPoints();
        }

        iteration++;

        if (moved_points < (data.n_points * MIN_PTS_MVMT_PCT))
        {
            convergence_criterion = true;
        }

        it_end = chrono::high_resolution_clock::now();

        elapsed = it_end - it_start;
        elapsed_time = elapsed.count();

        iter_t = elapsed_time;

        storeIterationData(file, data, clusters);


        // PRINTING
        #if PRINT >= VERBOSE
        cout << "ITERATION " << iteration-1 << endl;
        printClusters(clusters);
        printMovedPoints(moved_points, MIN_PTS_MVMT_PCT, iteration, convergence_criterion);
        cout << endl;
        #elif PRINT >= LOG
        printBenchmarkCSV(iteration-1, iter_t, dist_calc_t);
        #endif
    }

    // Trying out per-iteration data storage
    //storeData("clustered_data", data, clusters);

    // Benchmarking
    #if PRINT >= VERBOSE
    cout << "Total distance compute time: " << total_dist_calc_t << "(avg. per iteration: " << total_dist_calc_t/(iteration-1) << ")\n";
    #elif PRINT >= LOG
    //cout << "total dist. compute T" << ", " << "avg. T per iter" << endl;
    //cout << total_dist_calc_t << ", " << total_dist_calc_t/(iteration-1) << endl;
    #endif
}
