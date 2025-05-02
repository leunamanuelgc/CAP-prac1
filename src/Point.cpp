#include "Point.hpp"

Point::Point(int point_id, int cluster_id, int n_values, const std::vector<float> &values)
{
    this->point_id = point_id;
    this->cluster_id = cluster_id;
    this->values.clear();
    for (int i = 0; i < n_values; i++)
    {
        this->values.push_back(values[i]);
    }
}