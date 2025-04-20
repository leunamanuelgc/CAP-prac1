#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <random>

using namespace std;

#define PI 3.141582f

vector<float> getNBallPoint(vector<float> k, uint32_t dim, float max_radius, float min_radius = 0.0f) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float> norm_dist(0.0f, 1.0f);
    std::uniform_real_distribution<float> uniform_dist(0.0f, 1.0f);

    // 1. Generate a random direction vector (normal distribution)
    vector<float> direction(dim);
    for (int i = 0; i < dim; i++) {
        direction[i] = norm_dist(gen);
    }

    // 2. Normalize the direction vector
    float norm = std::sqrt(std::inner_product(direction.begin(), direction.end(), direction.begin(), 0.0f));
    for (float &val : direction) {
        val /= norm;
    }

    // 3. Generate radius r with uniform distribution in volume
    float u = uniform_dist(gen); // U[0,1]
    float r = dim < 2 ? u : std::pow(u, 1.0f / (dim-1.0)); // uniform in volume
    r = min_radius + r * (max_radius - min_radius);


    //cout << endl;
    // 4. Scale direction vector and add to center point
    vector<float> point(dim);
    for (int i = 0; i < dim; i++) {
        point[i] = k[i] + r * direction[i];
        //cout << point[i] << "\t";
    }

    return point;
}

int main()
{
    const int N_CLUSTERS = 10;
    const int N_POINTS_PER_CLUSTER = 500;
    const int N_DIMS = 3;

    const float MIN_VAL = 0.0;
    const float MAX_VAL = 35.0;

    FILE* resultsFile;
    resultsFile = fopen("salida", "wb");
    int nFilas = N_CLUSTERS * N_POINTS_PER_CLUSTER;
    int nCol = N_DIMS;
    fwrite(&nFilas, sizeof(int), 1, resultsFile);
    fwrite(&nCol, sizeof(int), 1, resultsFile);
    for (int i = 0; i < N_CLUSTERS; i++)
    {
        vector<float> centroid = getNBallPoint(vector<float>(N_DIMS, 0), N_DIMS, 100, 0);
        for (int j = 0; j < N_POINTS_PER_CLUSTER; j++)
            fwrite(getNBallPoint(centroid, N_DIMS, MAX_VAL, MIN_VAL).data(), sizeof(float) * nCol, 1, resultsFile);
    }

    fclose(resultsFile);
}