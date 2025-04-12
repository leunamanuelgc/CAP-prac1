#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <array>

using namespace std;

struct dataResult {
    int nFilas;
    int nColumnas;
    vector<float*> points;
};

dataResult readData(){
    FILE* inFile;
    dataResult data;
    int nFilas, nCol;
    //Abrir archivo para lectura en binario
    inFile = fopen("salida", "rb");
    if (inFile == NULL) {
        fputs("File error", stderr);
        exit(1);
    }
    //Leer el número de filas y columnas
    fread(&nFilas, sizeof(int),1, inFile);
    fread(&nCol, sizeof(int), 1, inFile);
    data.nFilas = nFilas;
    data.nColumnas = nCol;
    //Leer y guardar los puntos del fichero binario
    for (int i = 0; i<nFilas; i++){
        size_t pointSize = sizeof(float)*nCol;
        float* point = (float*) malloc(pointSize);
        size_t readCount = fread(point, sizeof(float), nCol, inFile);
        data.points.push_back(point);
    }
    //Cerrar archivo
    fclose(inFile);
    return data;
}

int main(int argc, char** argv){
    int K = 4;
    
    //Obtencion de los puntos
    dataResult data = readData();

    for (int i = 0; i < data.points.size(); i++){
        for (int j=0; j<data.nColumnas; j++){
            cout << data.points[i][j]<< "\t";
        }
        cout << "\n";
    }

    // //Algoritmo k-medias
    // for (int i = 0; i<dataPoints.size(); i+=dataPoints.size()/4){
    //     //
    // }
    return 0;
}