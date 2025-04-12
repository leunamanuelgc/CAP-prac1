#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

struct point2D {
    float x;
    float y;
};

vector<point2D> readData(){
    FILE* inFile;
    vector<point2D>data;
    int nFilas, nCol;

    inFile = fopen("salida", "rb");
    if (inFile == NULL) {
        fputs("File error", stderr);
        exit(1);
    }
    
    fread(&nFilas, sizeof(int),1, inFile);
    fread(&nCol, sizeof(int), 1, inFile);
    for (int i = 0; i<nFilas; i++){
        point2D point;
        fread(&point, sizeof(point), 1, inFile);
        data.push_back(point);
    }

    fclose(inFile);
    return data;
}

int main(int argc, char** argv){
    vector<point2D> dataPoints = readData();

    for (int i = 0; i < dataPoints.size(); i++)
        cout << dataPoints[i].x << "\t" << dataPoints[i].y << "\n";

    return 0;
}