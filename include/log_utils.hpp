#pragma once

#include <vector>
#include <iostream>
#include <iomanip>

inline void printLocalPointInfo(std::vector<PointRef> points, int iter_n, int n_ranks, int rank)
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
        PointRef& p = points[i];//(start_idx+i*points.size()/MAX_POINTS_PRINT)%points.size()];
        
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

inline void printMovedPoints(int n_moved, float convrg_pct, int iter, bool converged)
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

inline void printVector(double* v, size_t dim)
{
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "[";
    for (int i = 0; i < dim; i++)
    {
        std::cout << v[i];
        if(i<(dim-1))
            std::cout << " ";
        else
            std::cout << "]";
    }
}

inline void printVector(float* v, size_t dim)
{
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "[";
    for (int i = 0; i < dim; i++)
    {
        std::cout << v[i];
        if(i<(dim-1))
            std::cout << " ";
        else
            std::cout << "]";
    }
}