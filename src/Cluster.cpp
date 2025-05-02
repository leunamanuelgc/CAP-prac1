#include "Cluster.hpp"

#include <algorithm>

Cluster::Cluster(int cluster_id, int n_dim)
{
    this->cluster_id = cluster_id;
    centroid = std::vector<float>(n_dim, 0);
}

Cluster::Cluster(int cluster_id, const std::vector<float> &centroid)
{
    this->cluster_id = cluster_id;
    this->centroid = centroid;
}

void Cluster::addPoint(const Point &np) { points.push_back(np); }

void Cluster::removePointAt(int index) { points.erase(points.begin() + index); }

void Cluster::removePointByID(int point_id)
{
    // in theory, the oldest poinst in a cluster are probably correctly categorized.
    // following this logic points that are to be removed are probably closer to the end.
    // probably not a worth while optimization but it SHOULD actually improve performance,
    // with high point counts and after the first few iterations of K-Means.
    std::vector<Point>::reverse_iterator rit_rem = find_if(points.rbegin(), points.rend(),
                                                        [point_id](const Point &p)
                                                        {
                                                            return point_id == p.getID();
                                                        });
    // check if it was actually found
    if (rit_rem != points.rend())
        points.erase(next(rit_rem).base());
}