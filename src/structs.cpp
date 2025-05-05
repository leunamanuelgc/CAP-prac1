#include "structs.hpp"

#include <mpi.h>
#include <cstddef>
#include <string.h>

void centroid_diff_sum_function(void *in_bytes, void *inout_bytes, int *len, MPI_Datatype *datatype)
{
    const uint32_t dim = (*len - 2 * sizeof(uint32_t)) / (2 * sizeof(float));

    #pragma omp parallel
    {
        int i = 0;
        #pragma omp for nowait
        for (int i = 0; i < dim; i++)
        {
            float* in_ptr = &((float*)in_bytes)[i*sizeof(float)];
            float* out_ptr = &((float*)inout_bytes)[i*sizeof(float)];
            *out_ptr += *in_ptr;
        }

        #pragma single nowait
        {
            uint32_t* in_ptr = &((uint32_t*)in_bytes)[i*sizeof(float)];
            uint32_t* out_ptr = &((uint32_t*)inout_bytes)[i*sizeof(float)];
            *out_ptr += *in_ptr;
        }
        #pragma single nowait
        {
            uint32_t* in_ptr = &((uint32_t*)in_bytes)[i*sizeof(float)+sizeof(uint32_t)];
            uint32_t* out_ptr = &((uint32_t*)inout_bytes)[i*sizeof(float)+sizeof(uint32_t)];
            *out_ptr += *in_ptr;
        }
    }
}

void CentroidDiff::copyBytesIntoFlatBuff(std::byte* data_ptr) const
{
    size_t dim = add_points_sum.size();

    memcpy(data_ptr, add_points_sum.data(), dim * sizeof(float));
    data_ptr += dim * sizeof(float);
    memcpy(data_ptr, rem_points_sum.data(), dim * sizeof(float));
    data_ptr += dim * sizeof(float);
    memcpy(data_ptr, &add_points_count, sizeof(uint32_t));
    data_ptr += sizeof(uint32_t);
    memcpy(data_ptr, &rem_points_count, sizeof(uint32_t));
}

void CentroidDiff::copyFromFlatBytes(const std::byte* flat_byte_ptr)
{
    //const uint32_t dim = CentroidDiff::dim(len);
    const uint32_t dim = add_points_sum.size();
    
    //add_points_sum.resize(dim);
    //rem_points_sum.resize(dim);

    memcpy(add_points_sum.data(), flat_byte_ptr, dim * sizeof(float));
    flat_byte_ptr += dim * sizeof(float);
    memcpy(rem_points_sum.data(), flat_byte_ptr, dim * sizeof(float));
    flat_byte_ptr += dim * sizeof(float);
    add_points_count = *(uint32_t*)flat_byte_ptr;
    flat_byte_ptr += sizeof(uint32_t);
    rem_points_count = *(uint32_t*)flat_byte_ptr;
}
