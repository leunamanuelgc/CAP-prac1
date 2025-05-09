#pragma once

#include "v_ops.hpp"

#include <vector>
#include <string.h>

struct PointRef
{
friend struct PointData;
private:
    const float *values;
    int *cluster_id;
    uint32_t point_id;
    uint32_t dim;

    PointRef(float *coords, int *cluster_id, uint32_t id, uint32_t dim) :
        values(coords),
        cluster_id(cluster_id),
        point_id(id),
        dim(dim) {}
public:
    inline uint32_t getID() const { return point_id; }
    inline int getClusterID() const { return *cluster_id; }
    inline uint32_t getDim() const { return dim; }
    inline const float* getValues() const { return values; }
    inline float getValueAt(const int idx) const { return values[idx]; }

    inline void setClusterID(int new_cluster_id) { *cluster_id = new_cluster_id; }
};

inline double sqrDist(const PointRef &p1, const PointRef &p2)
{
    return sqrDist(p1.getValues(), p2.getValues(), p1.getDim());
};