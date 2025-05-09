#include "PointData.hpp"
#include "CentroidDiffRef.hpp"

#include <mpi.h>
#include <cstddef>
#include <string.h>
#include <iostream>

void centroid_diff_sum_function(void *in_bytes, void *inout_bytes, int *len, MPI_Datatype *datatype)
{
    int type_bytes = 0;
    MPI_Type_size(*datatype, &type_bytes);
    const uint32_t dim = CentroidDiffRef::dim(type_bytes);

    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int i = 0; i < *len; i++)
        {
            // Add two flat_centroid_diffs
            CentroidDiffRef::addFlatDiffs((std::byte*)in_bytes + i*type_bytes, (std::byte*)inout_bytes + i*type_bytes, dim);
        }
    }
}

void CentroidDiffRef::addFlatDiffs(void *a, void *b, int dim)
{
    // ADD both add and rem sums
    for (int i = 0; i < dim*2; i++)
    {
        float* in_float_ptr = &(static_cast<float*>(a))[i];
        float* out_float_ptr = &(static_cast<float*>(b))[i];
        *out_float_ptr += *in_float_ptr;
    }

    uint32_t* in_addcount_ptr = (uint32_t*)&(static_cast<float*>(a))[dim*2];
    uint32_t* out_addcount_ptr = (uint32_t*)&(static_cast<float*>(b))[dim*2];
    *out_addcount_ptr += *in_addcount_ptr;

    uint32_t* in_remcount_ptr = &((uint32_t*)&(static_cast<float*>(a))[dim*2])[1];
    uint32_t* out_remcount_ptr = &((uint32_t*)&(static_cast<float*>(b))[dim*2])[1];
    *out_remcount_ptr += *in_remcount_ptr;
}

void CentroidDiffRef::copyBytesIntoFlatBuff(std::byte* data_ptr) const
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

void CentroidDiffRef::copyFromFlatBytes(const std::byte* flat_byte_ptr)
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
