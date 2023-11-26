#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>

// 欧拉-黎曼Zeta函数计算函数
double euler_raman_zeta(int k, double s) {
    double total_sum = 0.0;

    #pragma omp parallel for reduction(+:total_sum)
    for (int i = 1; i <= k; i++) {
        for (int j = 1; j <= k; j++) {
            total_sum += pow(2, s) * pow(-1, i + j) / pow(i + j, s);
        }
    }

    return total_sum;
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <MAX_THREADS> <s> <k>\n", argv[0]);
        return 1;
    }

    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);

    int MAX_THREADS = atoi(argv[1]);
    double s = atof(argv[2]);
    int k = atoi(argv[3]);

    if (s <= 1.0 || k <= 0) {
        printf("s:%f must be greater than 1.0 and k must be a positive integer.\n", s);
        return 1;
    }

    int num_threads = (k < MAX_THREADS) ? k : MAX_THREADS;

    omp_set_num_threads(num_threads);

    double total_sum = euler_raman_zeta(k, s);

    gettimeofday(&end_time, NULL);
    double execution_time = (end_time.tv_sec - start_time.tv_sec) +
                           (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

    printf("Number of threads used: %d\n", num_threads);
    printf("Estimated value of ζ(%.2f) with k = %d is: %f\n", s, k, total_sum);
    printf("Execution time: %.6f seconds\n", execution_time);

    return 0;
}
