#include "structs.hpp"

#include <mpi.h>
#include <cstddef>
#include <string.h>

void centroid_diff_sum_function(void *in_bytes, void *inout_bytes, int *len, MPI_Datatype *datatype)
{
    const uint32_t dim = CentroidDiff::dim(*len);

    #pragma omp parallel
    {
        #pragma omp for nowait
        // ADD both add and rem sums
        for (int i = 0; i < dim*2; i++)
        {
            float* in_ptr = &(static_cast<float*>(in_bytes))[i];
            float* out_ptr = &(static_cast<float*>(inout_bytes))[i];
            *out_ptr += *in_ptr;
        }

        #pragma single nowait
        {
            uint32_t* in_ptr = (uint32_t*)&(static_cast<float*>(in_bytes))[dim*2];
            uint32_t* out_ptr = (uint32_t*)&(static_cast<float*>(inout_bytes))[dim*2];
            *out_ptr += *in_ptr;
        }
        #pragma single nowait
        {
            uint32_t* in_ptr = &((uint32_t*)&(static_cast<float*>(in_bytes))[dim*2])[1];
            uint32_t* out_ptr = &((uint32_t*)&(static_cast<float*>(inout_bytes))[dim*2])[1];
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
