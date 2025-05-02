#pragma once

#include "Point.hpp"

#include <vector>

#define NONE_CLUSTER -1

class Cluster
{
private:
    int cluster_id;
    std::vector<float> centroid;
    std::vector<Point> points;

public:
    //inline Cluster() {}
    inline Cluster(int n_dim) : Cluster(NONE_CLUSTER, n_dim) {}
    Cluster(int cluster_id, int n_dim);
    Cluster(int cluster_id, const std::vector<float> &centroid);

    inline int getID() const { return cluster_id; }
    inline int getDim() const { return centroid.size(); }
    inline std::vector<float> getCentroid() const { return centroid; }
    inline float getCentroidValue(int i) const { return centroid[i]; }
    inline int getNumPoints() const { return points.size(); }
    inline std::vector<Point> getPoints() const { return points; }
    inline Point getPointAt(int index) const { return points[index]; }

    inline void setCentroid(const std::vector<float> &new_centroid) { centroid = new_centroid; }
    inline void setCentroidValue(float value, int i) { centroid[i] = value; }

    // functions by copy to keep the point std::vector cache friendly, I think this is probably
    // for the best since we're going to be dividing clusters even further in memory with MPI.
    void addPoint(const Point &np);
    void removePointAt(int index);
    void removePointByID(int point_id);
};