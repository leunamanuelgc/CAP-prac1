#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

using namespace std;


vector<float> getRandomPoint(uint32_t dim, float min, float max)
{
    vector<float> p(dim);
    //cout << endl;
    for (int i = 0; i < dim; i++)
    {
        p[i] = min + (max-min) * (float)rand()/RAND_MAX;
        //cout << p[i] << "\t";
    }
    return p;
};

int main()
{
    const int N_POINTS = 10000;
    const int N_DIMS = 3;

    const float MIN_VAL = -100.0;
    const float MAX_VAL = 100.0;

    // seed rng 'randomly' to get different results per run
    srand(time(NULL));

    FILE* resultsFile;
    resultsFile = fopen("salida", "wb");
    int nFilas = N_POINTS;
    int nCol = N_DIMS;
    fwrite(&nFilas, sizeof(int), 1, resultsFile);
    fwrite(&nCol, sizeof(int), 1, resultsFile);

    for(int i = 0; i < N_POINTS; i++)
    {
        fwrite(getRandomPoint(N_DIMS, MIN_VAL, MAX_VAL).data(), sizeof(float) * nCol, 1, resultsFile);
    }
    
    fclose(resultsFile);
}