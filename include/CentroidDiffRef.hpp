#pragma once

struct CentroidDiffRef
{
private:
    friend struct CentroidDiffs;
    inline CentroidDiffRef(double* aps, double* rps, uint32_t* apc, uint32_t* rpc, uint32_t dim) :
        add_points_sum(aps),
        rem_points_sum(rps),
        add_points_count(apc),
        rem_points_count(rpc),
        dim(dim) {}
public:
    double *add_points_sum;
    double *rem_points_sum;
    uint32_t *add_points_count;
    uint32_t *rem_points_count;
    uint32_t dim;

    inline size_t nFlatBytes() const { return dim * sizeof(float) * 2 + sizeof(uint32_t) * 2; }
    inline static uint32_t dimsFromBytes(size_t n_flat_bytes) { return (n_flat_bytes - 2 * sizeof(uint32_t)) / (2 * sizeof(float)); }
    inline void incrementRemCount() { (*rem_points_count)++; }
    inline void incrementAddCount() { (*add_points_count)++; }

    // Function to flatten data into a byte buffer
    void copyBytesIntoFlatBuff(std::byte *data_ptr) const;
    void copyFromFlatBytes(const std::byte *flat_byte_ptr);
    static void addFlatDiffs(void *a, void *b, int dim);
};

// Sum operator used in the reduction of the centroid diffs
void centroid_diff_sum_function(void *inputBuffer, void *outputBuffer, int *len, MPI_Datatype *datatype);
