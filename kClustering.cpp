#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
//#include <omp.h>
#include <iomanip>
#include <chrono>

using namespace std;

#define NONE_CLUSTER -1

#define VERBOSE 2
#define LOG 1

#define PRINT LOG
//#define PRINT VERBOSE

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
    float getValueAt(const int idx) const { return values[idx]; }

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

void storeDataHeader(ofstream& file_stream, pointData data, vector<Cluster> clusters)
{
    uint32_t n_clusters = clusters.size();
    file_stream.write(reinterpret_cast<char*>(&n_clusters), sizeof(n_clusters));
    file_stream.write(reinterpret_cast<char*>(&data.n_points), sizeof(data.n_points));
    file_stream.write(reinterpret_cast<char*>(&data.n_dim), sizeof(data.n_dim));
}

void storeIterationData(ofstream& file_stream, pointData data, vector<Cluster> clusters)
{
    uint32_t n_clusters = clusters.size();
    // Store cluster centroid and id data
    for (int i = 0; i < n_clusters; i++)
    {
        file_stream.write(reinterpret_cast<char*>(&clusters[i].getCentroid()[0]), sizeof(float) * clusters[i].getDim());
        auto cluster_id = clusters[i].getID();
        file_stream.write(reinterpret_cast<char*>(&cluster_id), sizeof(clusters[i].getID()));
    }
    // Store point coordinates and cluster_id data
    for (int i = 0; i < data.n_points; i++)
    {
        file_stream.write(reinterpret_cast<char*>(&data.points[i].getValues()[0]), sizeof(float) * data.n_dim);
        auto cluster_id = data.points[i].getClusterID();
        file_stream.write(reinterpret_cast<char*>(&cluster_id), sizeof(cluster_id));
    }
}

void storeData(const string& filename, pointData data, vector<Cluster> clusters)
{
    ofstream file(filename, ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    storeDataHeader(file, data, clusters);
    storeIterationData(file, data, clusters);
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
    #pragma omp parallel for default(shared) reduction(+:sqr_dist)
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

void printClusters(vector<Cluster> clusters)
{
    const int MAX_PRINT_LEN = 5;
    for (int i = 0; i < clusters.size(); i++)
    {
        cout << setprecision(6);
        cout << "Cluster[" << i << "]\n\tCentroid:\t";
        for (int j = 0; j < clusters[i].getDim(); j++)
        {
            cout << clusters[i].getCentroid()[j] << "\t";
        }
        cout << setprecision(2);
        cout <<  "\n\tPoints(" << clusters[i].getNumPoints() << "):\n\t\t{";
        for (int j = 0; j < clusters[i].getNumPoints(); j++)
        {
            if (j >= MAX_PRINT_LEN && j < (clusters[i].getNumPoints()-MAX_PRINT_LEN))
            {
                if (j == MAX_PRINT_LEN) cout << "...\n\t\t... ";
                continue;
            }
            cout << "[";
            for (int k = 0; k < clusters[i].getDim(); k++)
            {
                cout << clusters[i].getPointAt(j).getValueAt(k);
                if (k < (clusters[i].getDim()-1)) cout << ",";
            }
            cout << "]";
            if (j < (clusters[i].getNumPoints()-1))
                cout << " ";
        }
        cout << "}\n\n";

    }
};

void printMovedPoints(int n_moved, int convrg_pct, int iter, bool converged)
{
    if (converged)
    {
        cout << "n Moved Points(=" << n_moved << ") < "<< convrg_pct*100 << "% - convergence criterion met STOPPING K-MEANS after \n\t"
        <<  iter-1 << " iterations" << endl;
    }
    else
    {
        cout << "n Moved Points: " << n_moved << endl;
    }
}

void kMeans(pointData data, int k)
{
    const int MAX_ITERATIONS = 30;
    const float MIN_PTS_MVMT_PCT = 0.05;
    // Benchmarking
	auto start = chrono::high_resolution_clock::now();
	auto end = start;
	chrono::duration<double> elapsed;
	double elapsed_time;
    double total_dist_calc_t;
    
    
    vector<Cluster> clusters(k);
    cout << fixed << setprecision(6);

    // Inicializar los clusters
    for (int i = 0; i < k; i++)
    {
        // Pick out equally spaced points from the data input vector (i.e random points, except not bc yes)
        int centroid_idx = (data.n_points / k) * i + data.n_points/k/2;

        // RANDOM INIT
        //int centroid_idx = i;

        // RANDOM INIT
        //int centroid_idx = (int)((float)rand()/RAND_MAX * data.n_points);\
        //cout << "AAAAAAAAAAAAAAAAA: \t" << centroid_idx << endl;
        clusters[i] = Cluster(i, data.points[centroid_idx].getValues());
    }

    // Open file for iteration output
    ofstream file("clustered_data", ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: clustered_data");
    }
    // And write header
    storeDataHeader(file, data, clusters);

    int iteration = 0;
    bool convergence_criterion = false;
    while (iteration < MAX_ITERATIONS && !convergence_criterion)
    {
        int moved_points = 0;

        double dist_calc_t  = 0;

        for (int i = 0; i < data.n_points; i++)
        {
            Point& p = data.points[i];
            double min_dist = numeric_limits<double>::infinity();
            int closest_cluster_id;

            start = chrono::high_resolution_clock::now();
            // call the algorithm #1 here
            // Calcular las distancias entre los puntos y los centroides de los k nodos
            // *emabarrassingly parallel
            #pragma omp parallel for
            for (int j = 0; j < k; j++)
            {
                double dist = sqrDist(p.getValues(), clusters[j].getCentroid());
                
                #pragma omp critical
                {
                    // Comprobar cual es la menor distancia entre un punto y los centroides
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        closest_cluster_id = j;
                    }
                }
            }
            end = chrono::high_resolution_clock::now();
            elapsed = end - start;
            elapsed_time = elapsed.count();
            dist_calc_t += elapsed_time;

            //  sync

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

        // sync -> reduction (average)

        // Readjust new cluster centroid (cluster's points' average)
        for (int i = 0; i < k; i++)
        {
            if (clusters[i].getNumPoints() == 0) continue;

            vector<float> new_centroid(clusters[i].getDim(), 0);

            for (int j = 0; j < clusters[i].getNumPoints(); j++)
            {
                for (int l = 0; l < new_centroid.size(); l++)
                {
                    new_centroid[l] += clusters[i].getPointAt(j).getValueAt(l);
                }
            }
            for (int j = 0; j < new_centroid.size(); j++)
            {
                new_centroid[j] /= clusters[i].getNumPoints();
            }
            clusters[i].setCentroid(new_centroid);
        }

        iteration++;

        if (moved_points < (data.n_points * MIN_PTS_MVMT_PCT))
        {
            convergence_criterion = true;
        }

        total_dist_calc_t += dist_calc_t;

        storeIterationData(file, data, clusters);


        // PRINTING
        #if PRINT >= VERBOSE
        cout << "iteration " << iteration << endl;
        printClusters(clusters);
        printMovedPoints(moved_points, MIN_PTS_MVMT_PCT, iteration, convergence_criterion);
        cout << endl;
        #elif PRINT >= LOG
        cout << "[" << iteration-1 <<"] dist. T, " << dist_calc_t << endl;
        #endif
    }

    // Trying out per-iteration data storage
    //storeData("clustered_data", data, clusters);

    // Benchmarking
    #if PRINT >= VERBOSE
    cout << "Total distance compute time: " << total_dist_calc_t << "(avg. per iteration: " << total_dist_calc_t/(iteration-1) << ")\n";
    #elif PRINT >= LOG
    //cout << "total dist. compute T" << ", " << "avg. T per iter" << endl;
    //cout << total_dist_calc_t << ", " << total_dist_calc_t/(iteration-1) << endl;
    #endif
}

bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

int main(int argc, char **argv)
{
    const int N_CLUSTERS = 8;

    // Obtencion de los puntos
    pointData data = readData("salida");

    // Hacer backup de los datos de salida preexistentes
    if (fileExists("clustered_data") && argv[0] == "b")
    {
        int backup_id = 0;
        bool file_created = false;
        while (!file_created)
        {
            string backup_name = "clustered_data_" + to_string(backup_id);
            if (!fileExists(backup_name))
            {
                file_created = (rename("clustered_data", backup_name.c_str()) == 0);
            }
            backup_id++;
        }
    }

    //srand(time(NULL));

    kMeans(data, N_CLUSTERS);

    return 0;
}