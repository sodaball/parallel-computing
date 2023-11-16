#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#define NUM_POINTS 150
#define NUM_FEATURES 4
#define NUM_CLUSTERS 3
#define MAX_ITERATIONS 100

double data[NUM_POINTS][NUM_FEATURES]; // 保存数据点
double centroids[NUM_CLUSTERS][NUM_FEATURES]; // 保存聚类中心
int cluster_assignment[NUM_POINTS]; // 保存数据点分配的聚类

// 初始化聚类中心
void initialize_centroids() {
    // 随机选择数据点作为初始聚类中心
    for (int i = 0; i < NUM_CLUSTERS; i++) {
        int rand_idx = rand() % NUM_POINTS;
        memcpy(centroids[i], data[rand_idx], sizeof(double) * NUM_FEATURES);
    }
}

// 计算两个数据点之间的欧氏距离的平方
double calculate_distance_square(double* point1, double* point2) {
    double distance = 0.0;
    for (int i = 0; i < NUM_FEATURES; i++) {
        distance += pow(point1[i] - point2[i], 2);
    }
    return distance;
}

// 分配数据点到最近的聚类中心
void assign_clusters() {
    #pragma omp parallel for
    for (int i = 0; i < NUM_POINTS; i++) {
        double min_distance = calculate_distance_square(data[i], centroids[0]);
        int cluster = 0;

        for (int j = 1; j < NUM_CLUSTERS; j++) {
            double distance = calculate_distance_square(data[i], centroids[j]);
            if (distance < min_distance) {
                min_distance = distance;
                cluster = j;
            }
        }

        cluster_assignment[i] = cluster;
    }
}

// 更新聚类中心
void update_centroids() {
    #pragma omp parallel for
    for (int i = 0; i < NUM_CLUSTERS; i++) {
        for (int j = 0; j < NUM_FEATURES; j++) {
            double sum = 0.0;
            int count = 0;

            for (int k = 0; k < NUM_POINTS; k++) {
                if (cluster_assignment[k] == i) {
                    sum += data[k][j];
                    count++;
                }
            }

            if (count > 0) {
                centroids[i][j] = sum / count;
            }
        }
    }
}

int main() {
    // 读取Iris数据集
    FILE *fp = fopen("./data/iris.data", "r");
    if (fp == NULL) {
        printf("Error: Could not open the data file.\n");
        return 1;
    }

    for (int i = 0; i < NUM_POINTS; i++) {
        for (int j = 0; j < NUM_FEATURES; j++) {
            if (fscanf(fp, "%lf,", &data[i][j]) == EOF) {
                printf("Error: Insufficient data in the data file.\n");
                fclose(fp);
                return 1;
            }
        }
    }

    fclose(fp);

    /* 并行计算 */
    int num_threads = 4;
    omp_set_num_threads(num_threads);

    initialize_centroids();

    double start_time = omp_get_wtime();

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        assign_clusters();
        update_centroids();
    }

    double end_time = omp_get_wtime();

    // // 打印每个数据点的聚类分配结果
    // for (int i = 0; i < NUM_POINTS; i++) {
    //     printf("Data point %d is in cluster %d\n", i, cluster_assignment[i]);
    // }

    printf("the number of thread : %d, execution time: %f seconds\n", num_threads, end_time - start_time);

    /* 并行计算 */
    num_threads = 2;
    omp_set_num_threads(num_threads);

    initialize_centroids();

    start_time = omp_get_wtime();

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        assign_clusters();
        update_centroids();
    }

    end_time = omp_get_wtime();

    // // 打印每个数据点的聚类分配结果
    // for (int i = 0; i < NUM_POINTS; i++) {
    //     printf("Data point %d is in cluster %d\n", i, cluster_assignment[i]);
    // }

    printf("the number of thread : %d, execution time: %f seconds\n", num_threads, end_time - start_time);

    /* 串行计算 */
    num_threads = 1;
    omp_set_num_threads(num_threads);

    initialize_centroids();

    start_time = omp_get_wtime();

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        assign_clusters();
        update_centroids();
    }

    end_time = omp_get_wtime();

    // // 打印每个数据点的聚类分配结果
    // for (int i = 0; i < NUM_POINTS; i++) {
    //     printf("Data point %d is in cluster %d\n", i, cluster_assignment[i]);
    // }

    printf("the number of thread : %d, execution time: %f seconds\n", num_threads, end_time - start_time);

    return 0;
}
