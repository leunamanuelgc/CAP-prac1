#pragma once

#include "Cluster.hpp"

#include <vector>
#include <iostream>
#include <iomanip>

inline void printClusters(std::vector<Cluster> clusters) 
{
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
        std::cout <<  "\n\tPoints(" << clusters[i].getNumPoints() << "):\n\t\t{";
        for (int j = 0; j < clusters[i].getNumPoints(); j++)
        {
            if (j >= MAX_PRINT_LEN && j < (clusters[i].getNumPoints()-MAX_PRINT_LEN))
            {
                if (j == MAX_PRINT_LEN) std::cout << "...\n\t\t... ";
                continue;
            }
            std::cout << "[";
            for (int k = 0; k < clusters[i].getDim(); k++)
            {
                std::cout << clusters[i].getPointAt(j).getValueAt(k);
                if (k < (clusters[i].getDim()-1)) std::cout << ",";
            }
            std::cout << "]";
            if (j < (clusters[i].getNumPoints()-1))
            std::cout << " ";
        }
        std::cout << "}\n\n";

    }
};

inline void printMovedPoints(int n_moved, int convrg_pct, int iter, bool converged) {
    if (converged)
    {
        std::cout << "n Moved Points(=" << n_moved << ") < "
        << convrg_pct*100 << "% - convergence criterion met STOPPING K-MEANS after \n\t"
        <<  iter-1 << " iterations" << std::endl;
    }
    else
    {
        std::cout << "n Moved Points: " << n_moved << std::endl;
    }
}

inline void printBenchmarkCSV(int iter_n, double iter_t, double dist_t) {
    std::cout
    << "[" << iter_n-1 <<"], "
    << "iter. T, " << iter_t << ", "
    << "dist. T, " << dist_t <<  ", "
    << std::endl;
}