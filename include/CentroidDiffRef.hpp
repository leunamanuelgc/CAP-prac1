#pragma once

struct CentroidDiffRef
{
private:
    friend struct CentroidDiffs;
    
    uint32_t *add_points_count;
    uint32_t *rem_points_count;

    inline CentroidDiffRef(double* aps, double* rps, uint32_t* apc, uint32_t* rpc, uint32_t dim) :
        add_points_sum(aps),
        rem_points_sum(rps),
        add_points_count(apc),
        rem_points_count(rpc),
        dim(dim) {}
public:
    double *add_points_sum;
    double *rem_points_sum;
    const uint32_t dim;


    inline auto getAddPointsCount() const { return *add_points_count; }
    inline auto getRemPointsCount() const { return *rem_points_count; }
    inline void incrementRemCount() { (*rem_points_count)++; }
    inline void incrementAddCount() { (*add_points_count)++; }
};
