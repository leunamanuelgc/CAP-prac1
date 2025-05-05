#pragma once

#include "v_ops.hpp"

#include <vector>
#include <string.h>

class Point
{
private:
    std::vector<float> values;
    uint32_t point_id;
    int cluster_id;

public:
    inline Point() {}
    inline Point(uint32_t dim) : values(dim, 0), point_id(-1), cluster_id(-1) {}
    inline Point(uint32_t point_id, int cluster_id, const std::vector<float> &values) : point_id(point_id), cluster_id(cluster_id), values(values) {}

    inline int getID() const { return point_id; }
    inline int getClusterID() const { return cluster_id; }
    inline int getDim() const { return values.size(); }
    inline std::vector<float> getValues() const { return values; }
    inline float getValueAt(const int idx) const { return values[idx]; }

    inline void setID(int new_id) { point_id = new_id; }
    inline void setClusterID(int new_cluster_id) { cluster_id = new_cluster_id; }
    void copyValueMemData(const std::vector<float> new_value_list, uint32_t first, uint32_t n_values);
};

inline double sqrDist(const Point &p1, const Point &p2)
{
    return sqrDist(p1.getValues(), p2.getValues());
};