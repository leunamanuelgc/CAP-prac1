#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

/*
Compile with:
g++ generateRndSeeds.cpp -o generate.exe
After that execute with:
./generate.exe
*/

using namespace std;

#define PI 3.141582f

struct point3D {
    float x;
    float y;
    float z;
};

point3D getRandomPoint(float x0, float y0, float z0, float maxRadius, float minRadius = 0.0f )
{
    point3D p;
    float r = minRadius + (maxRadius-minRadius) * (float)rand() / RAND_MAX;
    float theta = 2.0f*PI*(float)rand()/RAND_MAX;
    float phi = acos(2.0f * (float)rand()/RAND_MAX - 1.0f);
    p.x = x0 + r * sin(phi) * cos(theta);
    p.y = y0 + r * sin(phi) * sin(theta);
    p.z = z0 + r * cos(phi);
    return p;
};

int main()
{
    int nClusters = 1;
    int nPointsPerCluster = 50;

    std::vector<point3D> data;
    for (int i = 0; i < nClusters; i++)
    {
        point3D centroid = getRandomPoint(0.0f, 0.0f, 0.0f, 20.0, 0.0);
        for (int j = 0; j < nPointsPerCluster; j++)
            data.push_back(getRandomPoint(centroid.x,centroid.y,centroid.z, 1.0f));
    }
    FILE* resultsFile;
    resultsFile = fopen("salida", "wb");
    int nFilas = nClusters * nPointsPerCluster;
    int nCol = 3;
    fwrite(&nFilas, sizeof(int), 1, resultsFile);
    fwrite(&nCol, sizeof(int), 1, resultsFile);
    fwrite(data.data(), sizeof(float), data.size()*nCol, resultsFile);
    fclose(resultsFile);
    for (int i = 0; i < data.size(); i++)
        std::cout << data[i].x << "\t" << data[i].y << "\t" << data[i].z << "\n";
}