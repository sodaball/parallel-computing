#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

struct ThreadArgs {
    int thread_id;
    double s;  // 将s参数定义为double类型
    int k_start;
    int k_end;
    double partial_sum;
};

// 欧拉-黎曼Zeta函数计算函数
void *euler_raman_zeta_thread(void *arg) {
    struct ThreadArgs *thread_args = (struct ThreadArgs *)arg;
    double sum = 0.0;
    for (int i = thread_args->k_start; i <= thread_args->k_end; i++) {
        for (int j = 1; j <= thread_args->k_end; j++) {
            sum += pow(2, thread_args->s) * pow(-1, i + j) / pow(i + j, thread_args->s);
        }
    }
    thread_args->partial_sum = sum;
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <MAX_THREADS> <s> <k>\n", argv[0]);
        return 1;
    }

    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);

    int MAX_THREADS = atoi(argv[1]);
    double s = atof(argv[2]); // 将s参数解析为double类型
    int k = atoi(argv[3]);

    if (s <= 1.0 || k <= 0) {
        printf("s:%f must be greater than 1.0 and k must be a positive integer.\n",s);
        return 1;
    }

    int num_threads = (k < MAX_THREADS) ? k : MAX_THREADS;
    pthread_t threads[num_threads];
    struct ThreadArgs thread_args[num_threads];
    double total_sum = 0.0;

    for (int i = 0; i < num_threads; i++) {
        thread_args[i].thread_id = i;
        thread_args[i].s = s; // 传递s参数给线程函数
        int chunk_size = k / num_threads;
        int start = i * chunk_size + 1;
        int end = (i == (num_threads - 1)) ? k : start + chunk_size - 1;
        thread_args[i].k_start = start;
        thread_args[i].k_end = end;

        int result = pthread_create(&threads[i], NULL, euler_raman_zeta_thread, &thread_args[i]);
        if (result) {
            printf("Error creating thread %d\n", i);
            exit(1);
        }
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
        total_sum += thread_args[i].partial_sum;
    }

    gettimeofday(&end_time, NULL);
    double execution_time = (end_time.tv_sec - start_time.tv_sec) +
                           (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

    printf("Number of threads used: %d\n", num_threads);
    printf("Estimated value of ζ(%.2f) with k = %d is: %f\n", s, k, total_sum);
    printf("Execution time: %.6f seconds\n", execution_time);

    return 0;
}
