#include "Point.hpp"

/// @brief READS RAW DATA (careful for pointer manipulation) Not very idiomatic but most efficient
///     to copy a contiguous range of float values (which follows the spec of vector<float>) into the same sized float vector
/// @param new_values_ptr
void Point::copyValueMemData(const std::vector<float> new_value_list, uint32_t first, uint32_t n_values)
{
    if (values.size() != n_values)
        throw std::invalid_argument("Not enough memory preallocated on Point[" + Point::point_id + (std::string)"] ");
    const float *new_values_ptr = &new_value_list[first];
    memcpy(values.data(), new_values_ptr, n_values * sizeof(float));
}
