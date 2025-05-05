#pragma once

#include "Cluster.hpp"

#include <vector>
#include <iostream>
#include <iomanip>

inline void printClusters(std::vector<Cluster> clusters)
{
    std::cout << std::fixed << std::setprecision(6);

    const int MAX_PRINT_LEN = 5;
    for (int i = 0; i < clusters.size(); i++)
    {
        std::cout << std::setprecision(6);
        std::cout << "Cluster[" << i << "]\n\tCentroid:\t";
        for (int j = 0; j < clusters[i].getDim(); j++)
        {
            std::cout << clusters[i].getCentroid()[j] << "\t";
        }
        std::cout << std::setprecision(2);
        std::cout << "\n\tPoints(" << clusters[i].getNumPoints() << "):\n\t\t{";
        for (int j = 0; j < clusters[i].getNumPoints(); j++)
        {
            if (j >= MAX_PRINT_LEN && j < (clusters[i].getNumPoints() - MAX_PRINT_LEN))
            {
                if (j == MAX_PRINT_LEN)
                    std::cout << "...\n\t\t... ";
                continue;
            }
            std::cout << "[";
            for (int k = 0; k < clusters[i].getDim(); k++)
            {
                std::cout << clusters[i].getPointAt(j).getValueAt(k);
                if (k < (clusters[i].getDim() - 1))
                    std::cout << ",";
            }
            std::cout << "]";
            if (j < (clusters[i].getNumPoints() - 1))
                std::cout << " ";
        }
        std::cout << "}\n\n";
    }
};

inline void printLocalPointInfo(std::vector<Point> points, int iter_n, int n_ranks, int rank)
{
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\nMPI_Rank[" << rank << "/" << n_ranks-1 << "] - n. Points: " << points.size()
    << "\t ITER(" << iter_n << ")" << std::endl;

    const int MAX_POINTS_PRINT = 10;
    int start_idx = rand()/RAND_MAX*points.size();

    std::cout << std::setprecision(2);
    auto print_p = std::min(points.size(), (size_t)MAX_POINTS_PRINT);
    for (int i = 0; i < print_p; i++)
    {
        Point& p = points[i];//(start_idx+i*points.size()/MAX_POINTS_PRINT)%points.size()];
        
        std::cout << "Point(" << p.getID() << ") => K(" << p.getClusterID() << "): \t";
        
        std::cout << "[";
        for (int j = 0; j < p.getDim(); j++)
        {
            std::cout << p.getValueAt(j);
            if (j < (p.getDim() - 1))
                std::cout << ",";
        }
        std::cout << "]" << std::endl;
    }
}

inline void printMovedPoints(int n_moved, int convrg_pct, int iter, bool converged)
{
    if (converged)
    {
        std::cout << "n Moved Points(=" << n_moved << ") < "
                  << convrg_pct * 100 << "% - convergence criterion met STOPPING K-MEANS after \n\t"
                  << iter - 1 << " iterations" << std::endl;
    }
    else
    {
        std::cout << "n Moved Points: " << n_moved << std::endl;
    }
}

inline void printBenchmarkCSV(int iter_n, double iter_t, double dist_t)
{
    std::cout
        << "[" << iter_n << "], "
        << "iter. T, " << iter_t << ", "
        << "dist. T, " << dist_t << std::endl;
}

inline void printVectors(const std::vector<float>& v, size_t dim)
{
    for (int i = 0; i < v.size(); i ++)
    {
        if (i % dim == 0) std::cout << "[";
        std::cout << (float)v[i];
        if ((i+1) % dim == 0) std::cout << "]" << std::endl;
        else std::cout << " ";
    }
}

inline void printVector(const std::vector<float>& v)
{
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "[";
    for (int i = 0; i < v.size(); i++)
    {
        std::cout << v[i];
        if(i<(v.size()-1))
            std::cout << " ";
        else
            std::cout << "]";
    }
}