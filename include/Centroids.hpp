#pragma once

struct Centroids
{
    std::vector<float> coords;
    uint32_t dim;
    int n;

public:
    inline Centroids(uint32_t n, uint32_t dim) : coords(n*dim, 0), dim(dim), n(n) {}

    inline float *operator [](int i) {return &coords[i*dim];}
    
};