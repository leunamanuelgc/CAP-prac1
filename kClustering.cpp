#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <iomanip>

using namespace std;

#define MAX_ITERATIONS 2000

class Point
{
private:
    int point_id;
    int cluster_id;
    int n_values;
    vector<float> values;

public:
    Point(int point_id, int cluster_id, int n_values, const vector<float> &values)
    {
        this->point_id = point_id;
        this->cluster_id = cluster_id;
        this->n_values = n_values;
        this->values.clear();
        for (int i = 0; i < n_values; i++)
        {
            this->values.push_back(values[i]);
        }
    };
    int getID() const { return point_id; }
    int getClusterID() const { return cluster_id; }
    int getNumValues() const { return n_values; }
    vector<float> getValues() const { return values; }

    bool operator==(const Point &other) const
    {
        for (int i = 0; i < n_values; i++)
        {
            if (values[i] != other.values[i])
                return false;
        }
        return true;
    }
};

struct dataResult
{
    int nFilas;
    int nColumnas;
    float **distances;
    vector<Point> points;
};

class Cluster
{
private:
    int cluster_id;
    vector<float> centroid;
    int n_values;
    vector<Point> points;

public:
    Cluster(int cluster_id, int n_values)
    {
        this->cluster_id = cluster_id;
        this->n_values = n_values;
        centroid.clear();
    }

    Cluster(int cluster_id, const vector<float> &centroid)
    {
        this->cluster_id = cluster_id;
        n_values = centroid.size();

        this->centroid.clear();
        for (int i = 0; i < n_values; i++)
        {
            this->centroid.push_back(centroid[i]);
        }
    }

    int getID() const { return cluster_id; }
    int getNumValues() const { return n_values; }
    vector<float> getCentroid() const { return centroid; }
    float getCentroidValue(int i) const { return centroid[i]; }
    int getNumPoints() const { return points.size(); }
    vector<Point> getPoints() const { return points; }
    Point getPointAt(int index) const { return points[index]; }

    void setCentroid(const vector<float> &new_centroid) { centroid = new_centroid; }
    void setCentroidValue(float value, int i) { centroid[i] = value; }
    // functions by copy to keep the point vector cache friendly, I think this is probably
    // for the best since we're going to be dividing clusters even further in memory with MPI.
    void addPoint(const Point &np) { points.push_back(np); }
    void removePointAt(int index) { points.erase(points.begin() + index); }
    void removePointByID(int point_id)
    {
        // in theory, the oldest poinst in a cluster are probably correctly categorized.
        // following this logic points that are to be removed are probably closer to the end.
        // probably not a worth while optimization but it SHOULD actually improve performance,
        // with high point counts and after the first few iterations of K-Means.
        auto rit_rem = find_if(points.rbegin(), points.rend(),
                               [point_id](const Point &p)
                               {
                                   return point_id == p.getID();
                               });
        points.erase(next(rit_rem).base());
    }
};

/// @brief Encargado de gestionar los clusters correspondientes a cada Worker-Node,
/// ademas del movimiento de puntos entre los mismos.
class Master_Node
{
};

/// @brief Encargado de ejecutar las tareas relevantes a los clusters que le corresponden.
class Worker_Node
{
private:
    //vector<Cluster> // Leaving this out for further versions, implementing only one cluster per node to start off.
    Cluster cluster;
    
public:
    Worker_Node(int data_dimensions, int data_points, vector<Point> data, int worker_id, )
};

/**
 * Read data from the file: "salida" and stores it inside a dataResult struct.
 * @returns dataResult
 */
dataResult readData()
{
    FILE *inFile;
    dataResult data;

    int nFilas, nCol;
    // Abrir archivo para lectura en binario
    inFile = fopen("salida", "rb");
    if (inFile == NULL)
    {
        fputs("File error", stderr);
        exit(1);
    }
    // Leer el número de filas y columnas
    fread(&nFilas, sizeof(int), 1, inFile);
    fread(&nCol, sizeof(int), 1, inFile);
    data.nFilas = nFilas;
    data.nColumnas = nCol;

    // Leer y guardar los puntos del fichero binario
    for (int i = 0; i < nFilas; i++)
    {
        
        std::vector<char> values(nCol, 0);
        std::ifstream input("salina.bin", std::ios::in | std::ifstream::binary);
        input.read(&values[0], values.size());

        Point point = new Point(i, clusterid, nCol, values);
        data.points.push_back(point);
    }
    // Cerrar archivo
    fclose(inFile);
    return data;
}

point initialize_centroid(int nColumnas)
{
    // Inicializar k centroides (aleatorios)
    point centroid;

    for (int i = 0; i < nColumnas; i++)
    {
        centroid.components.push_back(20.0f * rand() / RAND_MAX);
        cout << centroid.components[i] << "\t";
    }
    cout << "\n";
    return centroid;
}

/**
 * Calculates the squared euclidean distance of nDimensions, between the centroid and a point.
 * @param point
 * @param centroid
 * @param nColumnas Indicates the dimension of the points.
 * @returns Squared euclidean distance.
 */
float squared_euclidean_distance(point p, point centroid, int nColumnas)
{
    float distance = 0;
    for (int i = 0; i < nColumnas; i++)
    {
        distance += (p.components[i] - centroid.components[i]) * (p.components[i] - centroid.components[i]);
    }
    return distance;
}

/**
 * Returns the new centroid by calculating the mean of every component of the points
 * @param points Every point used to calculate the new centroid.
 * @param nColumnas Dimension of the points.
 * @returns New centroid.
 */
point calculate_new_centroid(vector<point> points, int nColumnas)
{
    point new_centroid;

    for (int i = 0; i < nColumnas; i++)
    {
        new_centroid.components.push_back(0.0);
        for (int j = 0; j < points.size(); j++)
        {
            new_centroid.components[i] += points[j].components[i];
        }
        new_centroid.components[i] /= points.size();
    }
    return new_centroid;
}

bool equivalentPoints(point p1, point p2, int nColumnas)
{
    for (int i = 0; i < nColumnas; i++)
    {
        if (p1.components[i] != p2.components[i])
            return false;
    }
    return true;
}

point addDistance(point p, float d)
{
    p.distance = d;
    return p;
}

void kMeans(dataResult data, int k)
{
    Node nodes[k];
    cout << fixed << setprecision(6);

    // Inicializar los centroides
    for (int i = 0; i < k; i++)
    {
        cout << "\t" << "centroid " << i << endl;
        nodes[i] = Node(initialize_centroid(data.nColumnas));
    }

    for (int iter = 0; iter < MAX_ITERATIONS; iter++)
    {
        cout << "iteration " << iter << endl;

        for (int i = 0; i < data.points.size(); i++)
        {
            float *sqr_distances = (float *)malloc(sizeof(float) * k);

            // Calcular las distancias entre los puntos y los centroides de los k nodos
            for (int j = 0; j < k; j++)
            {
                sqr_distances[j] = squared_euclidean_distance(data.points[i], nodes[j].centroid, data.nColumnas);
            }

            // Comprobar cual es la menor distancia entre un punto y los centroides
            float min_dist = numeric_limits<float>::infinity();
            int min_id = -1;
            for (int j = 0; j < k; j++)
            {
                if (sqr_distances[j] < min_dist)
                {
                    min_dist = sqr_distances[j];
                    min_id = j;
                }
            }

            // Buscar el punto dentro de los otros nodos. Si se encuentra, se elimina de ese nodo
            for (int j = 0; j < k; j++)
            {
                vector<point>::iterator it;
                // Si el vector de puntos no está vacío...

                if (!nodes[j].points.empty())
                {
                    it = nodes[j].findPoint(data.points[i]);

                    int last = nodes[j].points.size() - 1;
                    if (min_id == j)
                    {
                        if (it != nodes[j].points.end() || equivalentPoints(data.points[i], nodes[j].points[last], data.nColumnas))
                        {
                            // Si encuentra el punto en el mismo nodo en el que estaba, actualiza la distancia
                            int p_id = distance(nodes[j].points.begin(), it);
                            nodes[j].setDistanceInPoint(p_id, min_dist);
                        }
                        else
                        {
                            // Si no encuentra el punto en el nodo cuya distancia entre estos es menor, lo añade
                            addDistance(data.points[i], min_dist);
                            nodes[min_id].addPoint(data.points[i]);
                        }
                    }
                    else
                    {
                        if (it != nodes[j].points.end() || equivalentPoints(data.points[i], nodes[j].points[last], data.nColumnas))
                        {
                            // Si encuentra el punto en el nodo que no le pertenece, lo elimina
                            nodes[j].removePointByIter(it);
                            break;
                        }
                    }
                    // Si el vector de puntos está vacío...
                }
                else
                {
                    if (min_id == j)
                    {
                        nodes[min_id].addPoint(data.points[i]);
                        it = nodes[min_id].points.end();
                        int p_id = distance(nodes[min_id].points.begin(), it) - 1;
                        nodes[min_id].setDistanceInPoint(p_id, min_dist);
                    }
                }
            }
        }

        for (int i = 0; i < k; i++)
        {
            if (!nodes[i].points.empty())
                nodes[i].setCentroid(calculate_new_centroid(nodes[i].points, data.nColumnas));
            cout << "centroid " << i << "\t";
            for (int j = 0; j < data.nColumnas; j++)
            {
                cout << nodes[i].centroid.components[j] << "\t";
            }
            cout << "\n";
        }
        cout << endl;
    }

    // Asignar punto a su centroide mas cercano. Para ello, calcular la distancia entre un punto y todos los nodos.

    // Despues, cada nodo ajusta la media de todos los puntos

    // Repetir proceso hasta que el total de puntos desplazados sea menor del 5% o se hayan alcanzado 2000 iteraciones
}

int main(int argc, char **argv)
{
    int k = 4;

    // Obtencion de los puntos
    dataResult data = readData();

    // Reserva el espacio que usarán las distancais de los puntos
    data.distances = (float **)malloc(sizeof(float *) * k);
    for (int i = 0; i < k; i++)
    {
        data.distances[i] = (float *)malloc(sizeof(float) * data.nFilas);
    }

    // for (int i = 0; i < data.points.size(); i++){
    //     for (int j=0; j<data.nColumnas; j++){
    //         cout << data.points[i].components[j] << "\t";
    //     }
    //     cout << "\n";
    // }

    kMeans(data, k);

    return 0;
}