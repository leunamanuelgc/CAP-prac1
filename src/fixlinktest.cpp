#include <cstdint>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>


/**
 * Calculates the squared euclidean distance of N Dimensions, between two points.
 * @param p1
 * @param p2
 * @returns Squared euclidean distance.
 */
inline double sqrDist(const std::vector<float>& v1, const std::vector<float>& v2) {
    auto n_dim = v1.size();
    if (n_dim != v2.size())
        throw std::invalid_argument("Can't measure distance between points with different dimensions [P1-%d <=> P2-%d]");

    double sqr_dist = 0;
    #pragma omp parallel for default(shared) reduction(+:sqr_dist)
    for (int i = 0; i < n_dim; i++)
    {
        sqr_dist += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    
    return sqr_dist;
};

class Point
{
private:
    int point_id;
    int cluster_id;
    std::vector<float> values;

public:
    inline Point() {}
    Point(int point_id, int cluster_id, int n_values, const std::vector<float> &values);
    
    inline int getID() const { return point_id; }
    inline int getClusterID() const { return cluster_id; }
    inline int getDim() const { return values.size(); }
    inline std::vector<float> getValues() const { return values; }
    inline float getValueAt(const int idx) const { return values[idx]; }

    inline void setClusterID(int new_cluster_id) { cluster_id = new_cluster_id; }
};

inline double sqrDist(const Point& p1, const Point& p2)
{    
    return sqrDist(p1.getValues(), p2.getValues());
};

Point::Point(int point_id, int cluster_id, int n_values, const std::vector<float> &values)
{
    this->point_id = point_id;
    this->cluster_id = cluster_id;
    this->values.clear();
    for (int i = 0; i < n_values; i++)
    {
        this->values.push_back(values[i]);
    }
}

#define NONE_CLUSTER -1

class Cluster
{
private:
    int cluster_id;
    std::vector<float> centroid;
    std::vector<Point> points;

public:
    //inline Cluster() {}
    inline Cluster(int n_dim) : Cluster(NONE_CLUSTER, n_dim) {}
    Cluster(int cluster_id, int n_dim);
    Cluster(int cluster_id, const std::vector<float> &centroid);

    inline int getID() const { return cluster_id; }
    inline int getDim() const { return centroid.size(); }
    inline std::vector<float> getCentroid() const { return centroid; }
    inline float getCentroidValue(int i) const { return centroid[i]; }
    inline int getNumPoints() const { return points.size(); }
    inline std::vector<Point> getPoints() const { return points; }
    inline Point getPointAt(int index) const { return points[index]; }

    inline void setCentroid(const std::vector<float> &new_centroid) { centroid = new_centroid; }
    inline void setCentroidValue(float value, int i) { centroid[i] = value; }

    // functions by copy to keep the point std::vector cache friendly, I think this is probably
    // for the best since we're going to be dividing clusters even further in memory with MPI.
    void addPoint(const Point &np);
    void removePointAt(int index);
    void removePointByID(int point_id);
};

Cluster::Cluster(int cluster_id, int n_dim)
{
    this->cluster_id = cluster_id;
    centroid = std::vector<float>(n_dim, 0);
}

Cluster::Cluster(int cluster_id, const std::vector<float> &centroid)
{
    this->cluster_id = cluster_id;
    this->centroid = centroid;
}

void Cluster::addPoint(const Point &np) { points.push_back(np); }

void Cluster::removePointAt(int index) { points.erase(points.begin() + index); }

void Cluster::removePointByID(int point_id)
{
    // in theory, the oldest poinst in a cluster are probably correctly categorized.
    // following this logic points that are to be removed are probably closer to the end.
    // probably not a worth while optimization but it SHOULD actually improve performance,
    // with high point counts and after the first few iterations of K-Means.
    std::vector<Point>::reverse_iterator rit_rem = find_if(points.rbegin(), points.rend(),
                                                        [point_id](const Point &p)
                                                        {
                                                            return point_id == p.getID();
                                                        });
    // check if it was actually found
    if (rit_rem != points.rend())
        points.erase(next(rit_rem).base());
}


struct PointData
{
    uint32_t n_points;
    uint32_t n_dim;
    std::vector<Point> points;
};

struct CentroidDiff
{
    std::vector<double> add_points_sum;
    int add_points_count;
    std::vector<double> rem_points_sum;
    int rem_points_count;

    //inline CentroidDiff() {}
    inline CentroidDiff(uint32_t dim) : add_points_count(0), rem_points_count(0)
    {
        add_points_sum = std::vector<double>(dim, 0);
        rem_points_sum = std::vector<double>(dim, 0);
    }
};

using namespace std;
/**
 * Read data from the file: "salida" and stores it inside a dataResult struct.
 * @returns dataResult
 */

PointData readData(const string& filename)
{
    ifstream file(filename, ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    PointData data;
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
};

void storeDataHeader(ofstream& file_stream, PointData data, vector<Cluster> clusters)
{
    uint32_t n_clusters = clusters.size();
    file_stream.write(reinterpret_cast<char*>(&n_clusters), sizeof(n_clusters));
    file_stream.write(reinterpret_cast<char*>(&data.n_points), sizeof(data.n_points));
    file_stream.write(reinterpret_cast<char*>(&data.n_dim), sizeof(data.n_dim));
}

void storeIterationData(ofstream& file_stream, PointData data, vector<Cluster> clusters)
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

void storeData(const string& filename, PointData data, vector<Cluster> clusters)
{
    ofstream file(filename, ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    storeDataHeader(file, data, clusters);
    storeIterationData(file, data, clusters);
}


inline bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

inline void printClusters(std::vector<Cluster> clusters) 
{
    const int MAX_PRINT_LEN = 5;
    for (int i = 0; i < clusters.size(); i++)
    {
        std::cout << std::setprecision(6);
        std::cout << "Cluster[" << i << "]\n\tCentroid:\t";
        for (int j = 0; j < clusters[i].getDim(); j++)
        {
            std::cout << clusters[i].getCentroid()[j] << "\t";
        }
        std::cout << std::setprecision(2);
        std::cout <<  "\n\tPoints(" << clusters[i].getNumPoints() << "):\n\t\t{";
        for (int j = 0; j < clusters[i].getNumPoints(); j++)
        {
            if (j >= MAX_PRINT_LEN && j < (clusters[i].getNumPoints()-MAX_PRINT_LEN))
            {
                if (j == MAX_PRINT_LEN) std::cout << "...\n\t\t... ";
                continue;
            }
            std::cout << "[";
            for (int k = 0; k < clusters[i].getDim(); k++)
            {
                std::cout << clusters[i].getPointAt(j).getValueAt(k);
                if (k < (clusters[i].getDim()-1)) std::cout << ",";
            }
            std::cout << "]";
            if (j < (clusters[i].getNumPoints()-1))
            std::cout << " ";
        }
        std::cout << "}\n\n";

    }
};

inline void printMovedPoints(int n_moved, int convrg_pct, int iter, bool converged) {
    if (converged)
    {
        std::cout << "n Moved Points(=" << n_moved << ") < "
        << convrg_pct*100 << "% - convergence criterion met STOPPING K-MEANS after \n\t"
        <<  iter-1 << " iterations" << std::endl;
    }
    else
    {
        std::cout << "n Moved Points: " << n_moved << std::endl;
    }
}

inline void printBenchmarkCSV(int iter_n, double iter_t, double dist_t) {
    std::cout
    << "[" << iter_n-1 <<"], "
    << "iter. T, " << iter_t << ", "
    << "dist. T, " << dist_t <<  ", "
    << std::endl;
}


void kMeans(PointData data, int k);

#define VERBOSE 2
#define LOG 1

#define PRINT LOG
//#define PRINT VERBOSE

int main(int argc, char **argv)
{
    std::vector<Cluster> foo(4, Cluster(10));
    for (auto elem : foo)
        std::cout << elem.getID() << ' ';

    const int N_CLUSTERS = 10;

    // Obtencion de los puntos
    PointData data = readData("./data/salida");

    // Hacer backup de los datos de salida preexistentes
    if (fileExists("./data/clustered_data") && argv[0] == "b")
    {
        int backup_id = 0;
        bool file_created = false;
        while (!file_created)
        {
            std::string backup_name = "./data/clustered_data_" + std::to_string(backup_id);
            if (!fileExists(backup_name))
            {
                file_created = (rename("./data/clustered_data", backup_name.c_str()) == 0);
            }
            backup_id++;
        }
    }

    srand(time(NULL));

    kMeans(data, N_CLUSTERS);

    return 0;
}

using namespace std;


void kMeans(PointData data, int k)
{
    const int MAX_ITERATIONS = 30;
    const float MIN_PTS_MVMT_PCT = 0.05;
    // Benchmarking
	auto start = chrono::high_resolution_clock::now();
	auto end = start;
	chrono::duration<double> elapsed;
	double elapsed_time;
    double total_dist_calc_t;
    double total_centroid_calc_t;
    double iter_t;
    
    std::vector<Cluster> clusters;//{k, Cluster(data.n_dim)};
    cout << fixed << setprecision(6);

    // Inicializar los clusters
    for (int i = 0; i < k; i++)
    {
        // Pick out equally spaced points from the data input vector (i.e random points, except not bc yes)
        int centroid_idx = (data.n_points / k) * i + data.n_points/k/2;

        // RANDOM INIT
        //int centroid_idx = i;

        // RANDOM INIT
        //int centroid_idx = (int)((float)rand()/RAND_MAX * data.n_points);
        //cout << "Centroid[" << i << "]: \t" << centroid_idx << endl;
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

    vector<uint32_t> old_cluster_p_count(k, 0);
    while (iteration < MAX_ITERATIONS && !convergence_criterion)
    {
        auto it_start = chrono::high_resolution_clock:: now();
        auto it_end = it_start;
        int moved_points = 0;

        double dist_calc_t  = 0;
        vector<CentroidDiff> centroidDiffs;//(k, CentroidDiff(data.n_dim));

        //#pragma omp parallel for default(private)
        for (int i = 0; i < data.n_points; i++)
        {
            Point& p = data.points[i];
            double min_dist = numeric_limits<double>::infinity();
            int closest_cluster_id;

            //start = chrono::high_resolution_clock::now();
            
            // Calcular las distancias entre los puntos y los centroides de los k nodos
            // *emabarrassingly parallel
            #pragma omp parallel shared(p,clusters, min_dist)
            {
                double th_min_dist = numeric_limits<double>::infinity();
                int th_closest_cluster_id;
            
                #pragma omp for nowait
                for (int j = 0; j < k; j++)
                {
                    double dist = sqrDist(p.getValues(), clusters[j].getCentroid());
                    if (dist < th_min_dist)
                    {
                        th_min_dist = dist;
                        th_closest_cluster_id = j;
                    }
                }

                // Bernat:
                // critical section extracted from thread group loop
                // to avoid reocurring locks between threads, since 
                // threads are handed out a portion of the loop to
                // run, the critical section would be hit constantly
                // greatly increasing wait times for the 'mutex' to clear
                #pragma omp critical
                {
                    // Comprobar cual es la menor distancia entre un punto y los centroides
                    if (th_min_dist < min_dist)
                    {
                        min_dist = th_min_dist;
                        closest_cluster_id = th_closest_cluster_id;
                    }
                }
            }
            /*end = chrono::high_resolution_clock::now();
            elapsed = end - start;
            elapsed_time = elapsed.count();
            dist_calc_t += elapsed_time;*/

            //  sync

            // Assign point to the new closest cluster
            if (p.getClusterID() != closest_cluster_id)
            {
                int old_cluster_id = p.getClusterID();
                // Remove from previous cluster
                if (old_cluster_id != NONE_CLUSTER)
                {
                    clusters[old_cluster_id].removePointByID(p.getID());
                    
                    // Centroid diff calculation
                    centroidDiffs[old_cluster_id].rem_points_count++;
                    for (int j = 0; j < centroidDiffs[old_cluster_id].rem_points_sum.size(); j++)
                    {
                        centroidDiffs[old_cluster_id].rem_points_sum[j] += p.getValueAt(j);
                    }
                }
                
                clusters[closest_cluster_id].addPoint(p);

                // Cendtroid diff calculation
                centroidDiffs[closest_cluster_id].add_points_count++;
                for (int j = 0; j < centroidDiffs[closest_cluster_id].add_points_sum.size(); j++)
                {
                    centroidDiffs[closest_cluster_id].add_points_sum[j] += p.getValueAt(j);
                }

                p.setClusterID(closest_cluster_id);
                moved_points++;
            }
        }

        // sync -> reduction (~average) (actually we will do a weithed
        //         mean calculation as a reduction method for the centroid)

        // Readjust new cluster centroid (cluster's points' average)
        #pragma omp parallel for default(private)
        for (int i = 0; i < k; i++)
        {
            if (centroidDiffs[i].add_points_count + centroidDiffs[i].rem_points_count == 0)
                continue;

            // copy old centroid
            vector<float> new_centroid(clusters[i].getCentroid());
            // maybe faster if initizlization is done inside the loop too
            for (int j = 0; j < new_centroid.size(); j++)
            {
                new_centroid[j] *= (float)old_cluster_p_count[i] / clusters[i].getNumPoints();
            }
            
            if (centroidDiffs[i].add_points_count > 0)
            {
                for (int j = 0; j < new_centroid.size(); j++)
                {
                    float add_points_mean_val = centroidDiffs[i].add_points_sum[j] / centroidDiffs[i].add_points_count;
                    new_centroid[j] +=
                        add_points_mean_val *
                        (float)centroidDiffs[i].add_points_count/clusters[i].getNumPoints();
                }
            }
            if (centroidDiffs[i].rem_points_count > 0)
            {
                for (int j = 0; j < new_centroid.size(); j++)
                {
                    float rem_points_mean_val = centroidDiffs[i].rem_points_sum[j] / centroidDiffs[i].rem_points_count;
                    new_centroid[j] -=
                        rem_points_mean_val *
                        (float)centroidDiffs[i].rem_points_count/old_cluster_p_count[i];
                }
            }
            
            clusters[i].setCentroid(new_centroid);
            old_cluster_p_count[i] = clusters[i].getNumPoints();
        }

        iteration++;

        if (moved_points < (data.n_points * MIN_PTS_MVMT_PCT))
        {
            convergence_criterion = true;
        }

        it_end = chrono::high_resolution_clock::now();

        elapsed = it_end - it_start;
        elapsed_time = elapsed.count();

        iter_t = elapsed_time;

        storeIterationData(file, data, clusters);


        // PRINTING
        #if PRINT >= VERBOSE
        cout << "ITERATION " << iteration-1 << endl;
        printClusters(clusters);
        printMovedPoints(moved_points, MIN_PTS_MVMT_PCT, iteration, convergence_criterion);
        cout << endl;
        #elif PRINT >= LOG
        printBenchmarkCSV(iteration-1, iter_t, dist_calc_t);
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
