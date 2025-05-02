#pragma once

#include "v_ops.hpp"

#include <vector>

class Point
{
private:
    int point_id;
    int cluster_id;
    std::vector<float> values;

public:
    inline Point() {}
    Point(int point_id, int cluster_id, int n_values, const std::vector<float> &values);
    
    inline int getID() const { return point_id; }
    inline int getClusterID() const { return cluster_id; }
    inline int getDim() const { return values.size(); }
    inline std::vector<float> getValues() const { return values; }
    inline float getValueAt(const int idx) const { return values[idx]; }

    inline void setClusterID(int new_cluster_id) { cluster_id = new_cluster_id; }
};


inline double sqrDist(const Point& p1, const Point& p2)
{    
    return sqrDist(p1.getValues(), p2.getValues());
};