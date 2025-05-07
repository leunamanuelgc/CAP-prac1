#include <iostream>
#include <vector>
#include <cstdio>

using namespace std;

vector<float> getDebugPoint(int id, int dim)
{
    vector<float> p(dim);
    for (int i = 0; i < dim; i++)
    {
        p[i] = id + (float)i/dim;
    }

    return p;
}

int main()
{
    const int N_POINTS = 50;
    const int N_DIMS = 3;

    FILE* resultsFile;
    resultsFile = fopen("./build/data/salida", "wb");
    int nFilas = N_POINTS;
    int nCol = N_DIMS;
    fwrite(&nFilas, sizeof(int), 1, resultsFile);
    fwrite(&nCol, sizeof(int), 1, resultsFile);

    for(int i = 0; i < N_POINTS; i++)
    {
        fwrite(getDebugPoint(i, nCol).data(), sizeof(float) * nCol, 1, resultsFile);
    }
    
    fclose(resultsFile);
}