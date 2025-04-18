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
#define NONE_CLUSTER -1

class Point
{
private:
    int point_id;
    int cluster_id;
    int n_dim;
    vector<float> values;

public:
    Point(int point_id, int cluster_id, int n_values, const vector<float> &values)
    {
        this->point_id = point_id;
        this->cluster_id = cluster_id;
        this->n_dim = n_values;
        this->values.clear();
        for (int i = 0; i < n_values; i++)
        {
            this->values.push_back(values[i]);
        }
    };
    int getID() const { return point_id; }
    int getClusterID() const { return cluster_id; }
    int getDim() const { return n_dim; }
    vector<float> getValues() const { return values; }

    void setClusterID(int new_cluster_id) { cluster_id = new_cluster_id; }
    bool operator==(const Point &other) const
    {
        for (int i = 0; i < n_dim; i++)
        {
            if (values[i] != other.values[i])
                return false;
        }
        return true;
    }
};

struct pointData
{
    uint32_t n_points;
    uint32_t n_dim;
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
    Cluster() {}
    Cluster(int cluster_id, int n_dim)
    {
        this->cluster_id = cluster_id;
        this->n_values = n_dim;
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
    int getDim() const { return n_values; }
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
        vector<Point>::reverse_iterator rit_rem = find_if(points.rbegin(), points.rend(),
                                                           [point_id](const Point &p)
                                                           {
                                                               return point_id == p.getID();
                                                           });
        // check if it was actually found
        if (rit_rem != points.rend())
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
    // vector<Cluster> // Leaving this out for further versions, implementing only one cluster per node to start off.
    Cluster cluster;

public:
    Worker_Node(int data_dimensions, int data_points, vector<Point> data, int worker_id);
};

/**
 * Read data from the file: "salida" and stores it inside a dataResult struct.
 * @returns dataResult
 */
pointData readData(const string& filename)
{
    ifstream file(filename, ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    pointData data;
    uint32_t n_rows, n_cols;
    
    // Leer el número de filas y columnas
    file.read(reinterpret_cast<char*>(&n_rows), sizeof(n_rows));
    file.read(reinterpret_cast<char*>(&n_cols), sizeof(n_cols));
    data.n_points = n_rows;
    data.n_dim = n_cols;

    // Leer y guardar los puntos del fichero binario
    for (int i = 0; i < n_rows; i++)
    {
        vector<float> values(n_cols);
        file.read(reinterpret_cast<char*>(&values[0]), n_cols * sizeof(float));
        Point p(i, NONE_CLUSTER, n_cols, values);
        data.points.push_back(move(p));
    }

    return data;
}

/*
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
*/

/**
 * Calculates the squared euclidean distance of nDimensions, between two points.
 * @param p1
 * @param p2
 * @returns Squared euclidean distance.
 */
double sqrDist(const vector<float>& v1, const vector<float>& v2) {
    int n_dim = v1.size();
    if (n_dim != v2.size())
        throw invalid_argument("Can't measure distance between points with different dimensions [P1-%d <=> P2-%d]");

    double sqr_dist = 0;
    for (int i = 0; i < n_dim; i++)
    {
        sqr_dist += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return sqr_dist;
};

double sqrDist(const Point& p1, const Point& p2) {    
    return sqrDist(p1.getValues(), p2.getValues());
};

/**
 * Returns the new centroid by calculating the mean of every component of the points
 * @param points Every point used to calculate the new centroid.
 * @param nColumnas Dimension of the points.
 * @returns New centroid.
 */
/*
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
*/

/*
bool equivalentPoints(point p1, point p2, int nColumnas)
{
    for (int i = 0; i < nColumnas; i++)
    {
        if (p1.components[i] != p2.components[i])
            return false;
    }
    return true;
}
*/

/*
point addDistance(point p, float d)
{
    p.distance = d;
    return p;
}
*/

void kMeans(pointData data, int k)
{
    vector<Cluster> clusters(k);
    cout << fixed << setprecision(6);

    // Inicializar los clusters
    for (int i = 0; i < k; i++)
    {
        cout << "\t" << "centroid " << i << endl;
        // Pick out equally spaced points from the data input vector (i.e random points)
        int centroid_idx = (data.n_points / k) * i + data.n_points/k/2;
        clusters[i] = Cluster(i, data.points[centroid_idx].getValues());
    }

    int iteration = 0;
    bool convergence_criterion = false;
    while (iteration < MAX_ITERATIONS && !convergence_criterion)
    {
        cout << "iteration " << iteration << endl;

        int moved_points = 0;

        // *emabarrassingly parallel
        for (int i = 0; i < data.n_points; i++)
        {
            Point& p = data.points[i];
            double min_dist = numeric_limits<double>::infinity();
            int closest_cluster_id;

            // Calcular las distancias entre los puntos y los centroides de los k nodos
            for (int j = 0; j < k; j++)
            {
                double dist = sqrDist(p.getValues(), clusters[j].getCentroid());
                
                // Comprobar cual es la menor distancia entre un punto y los centroides
                if (dist < min_dist)
                {
                    min_dist = dist;
                    closest_cluster_id = j;
                }
            }
            // BERNAT - no estoy seguro de si es esto lo que hay que hacer, hasta aqui he llegado.
            //          hay que implementar el algoritmo de k-medias y ya basicamente.
            //          No me complicaria con intentar separarlo en los worker_nodes
            //          de primeras, porque probablemente complique el codigo mas de lo
            //          necesario para nada.

            // Assign point to the new closest cluster
            if (p.getClusterID() != closest_cluster_id)
            {
                // Remove from previous cluster
                if (p.getClusterID() != NONE_CLUSTER)
                {
                    clusters[p.getClusterID()].removePointByID(p.getID());
                }

                p.setClusterID(closest_cluster_id);
                clusters[closest_cluster_id].addPoint(p);
                moved_points++;
            }
        }

        // Readjust new cluster centroid (cluster's points' average)
        for (int i = 0; i < k; i++)
        {
            vector<float> new_centroid(clusters[i].getDim());

            for (int j = 0; j < clusters[i].getNumPoints(); j++)
            {
                vector<float> point_values = clusters[i].getPointAt(j).getValues();
                for (int l = 0; l < new_centroid.size(); l++)
                {
                    new_centroid[l] += point_values[l];
                }
            }
            for (int j = 0; j < new_centroid.size(); j++)
            {
                new_centroid[j] /= clusters[i].getNumPoints();
            }
            clusters[i].setCentroid(new_centroid);

            cout << "Cluster[" << i << "] centroid:\t";
            for (int j = 0; j < clusters[i].getDim(); j++)
            {
                cout << clusters[i].getCentroid()[j] << "\t";
            }
            cout << "\n";
        }
        cout << endl;

        if (moved_points < data.n_points * 0.05)
        {
            convergence_criterion = true;
            cout << "n Moved Points(=" << moved_points << ") < 5% - convergence criterion met STOPPING K-MEANS" << endl;
        }
        else
        {
            cout << "n Moved Points: " << moved_points << endl;
        }

        iteration++;
    }
}

int main(int argc, char **argv)
{
    int k = 4;

    // Obtencion de los puntos
    pointData data = readData("salida");

    kMeans(data, k);

    return 0;
}