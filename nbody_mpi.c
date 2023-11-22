// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <mpi.h>

// #define G 6.67430e-11  // 万有引力常数
// #define tmax 1000       // 最大时间步数
// #define dt 1.0         // 时间步长
// #define convergence_threshold 1e-9  // 收敛阈值

// typedef struct {
//     double mass;
//     double x, y, z;
//     double vx, vy, vz;
// } Body;

// void initialize_bodies(Body* bodies, int num_bodies, int rank, int size) {
//     int body_per_proc = num_bodies / size;
//     int start = rank * body_per_proc;
//     int end = (rank == size - 1) ? num_bodies : start + body_per_proc;

//     // 设置相同的随机种子
//     srand(1234);  // 可以选择任何非负整数

//     for (int i = start; i < end; i++) {
//         bodies[i].mass = 1e24;  // 随机设置质量
//         bodies[i].x = rand() / (double)RAND_MAX;  // 随机初始位置
//         bodies[i].y = rand() / (double)RAND_MAX;
//         bodies[i].z = rand() / (double)RAND_MAX;
//         bodies[i].vx = 0.0;
//         bodies[i].vy = 0.0;
//         bodies[i].vz = 0.0;
//     }
// }

// void compute_force(Body* my_bodies, int my_start, int my_end, Body* all_bodies, int n) {
//     double softening = 1e-6;  // 软化参数

//     for (int i = my_start; i < my_end; i++) {
//         double fx = 0.0, fy = 0.0, fz = 0.0;

//         for (int j = 0; j < n; j++) {
//             if (i != j) {
//                 double dx = all_bodies[j].x - my_bodies[i - my_start].x;
//                 double dy = all_bodies[j].y - my_bodies[i - my_start].y;
//                 double dz = all_bodies[j].z - my_bodies[i - my_start].z;
//                 double r = sqrt(dx * dx + dy * dy + dz * dz + softening);  // 添加软化参数
//                 double force = (G * my_bodies[i - my_start].mass * all_bodies[j].mass) / (r * r * r);
//                 fx += force * dx;
//                 fy += force * dy;
//                 fz += force * dz;
//             }
//         }

//         my_bodies[i - my_start].vx += fx * dt / my_bodies[i - my_start].mass;
//         my_bodies[i - my_start].vy += fy * dt / my_bodies[i - my_start].mass;
//         my_bodies[i - my_start].vz += fz * dt / my_bodies[i - my_start].mass;
//     }
// }

// double compute_total_velocity(Body* my_bodies, int my_start, int my_end) {
//     double total_velocity = 0.0;

//     #pragma omp parallel for reduction(+:total_velocity)
//     for (int i = my_start; i < my_end; i++) {
//         total_velocity += sqrt(my_bodies[i - my_start].vx * my_bodies[i - my_start].vx +
//                                my_bodies[i - my_start].vy * my_bodies[i - my_start].vy +
//                                my_bodies[i - my_start].vz * my_bodies[i - my_start].vz);
//     }

//     return total_velocity;
// }

// int main(int argc, char** argv) {
//     MPI_Init(&argc, &argv);

//     int rank, size; // 进程编号和进程总数
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     int n = 100; // N体问题中的物体数量
//     int my_bodies_count = n / size; // 每个进程负责的物体数量

//     // 创建全局物体数组
//     Body* all_bodies = (Body*)malloc(sizeof(Body) * n);
    
//     /**
//      * @brief 并行版本
//      */
    
//     // 初始化所有物体的质量、位置和速度信息
//     initialize_bodies(all_bodies, n, rank, size);

//     // 创建每个进程负责的局部物体数组
//     Body* my_bodies = (Body*)malloc(sizeof(Body) * my_bodies_count);
//     // 从全局数组中复制每个进程负责的局部物体信息
//     int my_start = rank * my_bodies_count;
//     int my_end = my_start + my_bodies_count;
//     for (int i = my_start; i < my_end; i++) {
//         my_bodies[i - my_start] = all_bodies[i];
//     }

//     double total_velocity_change = 1.0;  // 初始值，保证至少一次循环

//     // 计时开始
//     double start_time, end_time;
//     start_time = MPI_Wtime();

//     for (int t = 0; t < tmax && total_velocity_change > convergence_threshold; t++) {
//         // 同步全局物体信息
//         MPI_Allgather(my_bodies, my_bodies_count, MPI_DOUBLE, all_bodies, my_bodies_count, MPI_DOUBLE, MPI_COMM_WORLD);

//         // 计算受力和更新位置
//         compute_force(my_bodies, 0, my_bodies_count, all_bodies, n);
//         for (int i = 0; i < my_bodies_count; i++) {
//             my_bodies[i].x += my_bodies[i].vx * dt;
//             my_bodies[i].y += my_bodies[i].vy * dt;
//             my_bodies[i].z += my_bodies[i].vz * dt;

//             // 更新位置时使用 all_bodies
//             my_bodies[i].vx = all_bodies[i].vx;
//             my_bodies[i].vy = all_bodies[i].vy;
//             my_bodies[i].vz = all_bodies[i].vz;
//         }

//         // 计算速度的变化并更新收敛值
//         double old_total_velocity_change = total_velocity_change;
//         total_velocity_change = compute_total_velocity(my_bodies, 0, n);

//         // // 输出每个时间步结束时的部分位置信息
//         // if (t % 1000 == 0) {
//         //     for (int i = 0; i < my_bodies_count; i++) {
//         //         printf("Rank %d - Body %d: x=%.4f, y=%.4f, z=%.4f\n", rank, i, my_bodies[i].x, my_bodies[i].y, my_bodies[i].z);
//         //     }
//         // }

//         // 输出收敛信息
//         printf("Convergence at time step %d: %.6f\n", t, total_velocity_change);

//         // 如果速度变化小于阈值，则认为系统已经收敛
//         if (fabs(old_total_velocity_change - total_velocity_change) < convergence_threshold) {
//             printf("Converged at time step %d\n", t);
//             break;
//         }
//     }

//     // 计时结束
//     end_time = MPI_Wtime();
//     double parallel_execution_time = end_time - start_time;

//     // 释放内存和MPI资源
//     free(all_bodies);
//     free(my_bodies);

//     // 输出并行执行时间
//     if (rank == 0) {
//         printf("Parallel execution time: %.4f seconds\n", parallel_execution_time);
//     }

//     /**
//      * @brief 串行版本
//      */

//     if (rank == 0) {
//         // 重新初始化全局物体数组
//         initialize_bodies(all_bodies, n, rank, size);

//         total_velocity_change = 1.0;  // 重置收敛值

//         // 计时开始
//         start_time = MPI_Wtime();

//         for (int t = 0; t < tmax && total_velocity_change > convergence_threshold; t++) {
//             // 计算受力和更新位置
//             compute_force(all_bodies, 0, n, all_bodies, n);
//             for (int i = 0; i < n; i++) {
//                 all_bodies[i].x += all_bodies[i].vx * dt;
//                 all_bodies[i].y += all_bodies[i].vy * dt;
//                 all_bodies[i].z += all_bodies[i].vz * dt;
//             }

//             // 计算速度的变化并更新收敛值
//             double old_total_velocity_change = total_velocity_change;
//             total_velocity_change = compute_total_velocity(my_bodies, 0, n);

//             // // 输出每个时间步结束时的部分位置信息
//             // if (t % 1000 == 0) {
//             //     for (int i = 0; i < my_bodies_count; i++) {
//             //         printf("Rank %d - Body %d: x=%.4f, y=%.4f, z=%.4f\n", rank, i, my_bodies[i].x, my_bodies[i].y, my_bodies[i].z);
//             //     }
//             // }

//             // 输出收敛信息
//             printf("Convergence at time step %d: %.6f\n", t, total_velocity_change);

//             // 如果速度变化小于阈值，则认为系统已经收敛
//             if (fabs(old_total_velocity_change - total_velocity_change) < convergence_threshold) {
//                 printf("Converged at time step %d\n", t);
//                 break;
//             }
//         }

//         // 计时结束
//         end_time = MPI_Wtime();
//         double serial_execution_time = end_time - start_time;

//         // 输出串行执行时间
//         printf("Serial execution time: %.4f seconds\n", serial_execution_time);

//         // 计算并输出加速比
//         double speedup = serial_execution_time / parallel_execution_time;
//         printf("Speedup: %.4f\n", speedup);
//     }

//     // 输出进程数、物体个数和迭代次数
//     if (rank == 0) {
//         printf("Number of processes: %d\n", size);
//         printf("Number of bodies: %d\n", n);
//         printf("Number of iterations: %d\n", tmax);
//     }

//     MPI_Finalize();

//     return 0;
// }


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define G 6.67430e-11  // 万有引力常数
#define tmax 1000       // 最大时间步数
#define dt 1.0          // 时间步长
#define convergence_threshold 1e-6  // 收敛阈值

typedef struct {
    double mass;
    double x, y, z;
    double vx, vy, vz;
} Body;

void initialize_bodies(Body* bodies, int num_bodies, int rank, int size) {
    int body_per_proc = num_bodies / size;
    int start = rank * body_per_proc;
    int end = (rank == size - 1) ? num_bodies : start + body_per_proc;

    srand(1234);  // 设置相同的随机种子

    for (int i = start; i < end; i++) {
        bodies[i].mass = 1e24;  // 随机设置质量
        bodies[i].x = rand() / (double)RAND_MAX;  // 随机初始位置
        bodies[i].y = rand() / (double)RAND_MAX;
        bodies[i].z = rand() / (double)RAND_MAX;
        bodies[i].vx = rand() / (double)RAND_MAX;
        bodies[i].vy = rand() / (double)RAND_MAX;
        bodies[i].vz = rand() / (double)RAND_MAX;
    }
}

void compute_force(Body* my_bodies, int my_start, int my_end, Body* all_bodies, int n) {
    double softening = 1e-6;  // 软化参数

    for (int i = my_start; i < my_end; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                double dx = all_bodies[j].x - my_bodies[i - my_start].x;
                double dy = all_bodies[j].y - my_bodies[i - my_start].y;
                double dz = all_bodies[j].z - my_bodies[i - my_start].z;
                double r = sqrt(dx * dx + dy * dy + dz * dz + softening);  // 添加软化参数
                double force = (G * my_bodies[i - my_start].mass * all_bodies[j].mass) / (r * r * r);
                fx += force * dx;
                fy += force * dy;
                fz += force * dz;
            }
        }

        my_bodies[i - my_start].vx += fx * dt / my_bodies[i - my_start].mass;
        my_bodies[i - my_start].vy += fy * dt / my_bodies[i - my_start].mass;
        my_bodies[i - my_start].vz += fz * dt / my_bodies[i - my_start].mass;
    }
}

double compute_total_velocity(Body* my_bodies, int my_start, int my_end) {
    double total_velocity = 0.0;

    #pragma omp parallel for reduction(+:total_velocity)
    for (int i = my_start; i < my_end; i++) {
        total_velocity += sqrt(my_bodies[i - my_start].vx * my_bodies[i - my_start].vx +
                               my_bodies[i - my_start].vy * my_bodies[i - my_start].vy +
                               my_bodies[i - my_start].vz * my_bodies[i - my_start].vz);
    }

    return total_velocity;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size; // 进程编号和进程总数
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 1000; // N体问题中的物体数量
    int my_bodies_count = n / size; // 每个进程负责的物体数量

    // 创建全局物体数组
    Body* all_bodies = (Body*)malloc(sizeof(Body) * n);

    // 初始化每个进程负责的局部物体数组
    Body* my_bodies = (Body*)malloc(sizeof(Body) * my_bodies_count);

    // 复制全局数组到每个进程的局部数组
    MPI_Scatter(all_bodies, my_bodies_count, MPI_DOUBLE, my_bodies, my_bodies_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double total_velocity_change = 1.0;  // 初始值，保证至少一次循环
    double start_time, end_time;

    start_time = MPI_Wtime();

    for (int t = 0; t < tmax && total_velocity_change > convergence_threshold; t++) {
        // 将每个进程的局部物体数组广播到所有进程
        MPI_Bcast(my_bodies, my_bodies_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // 计算受力和更新位置
        compute_force(my_bodies, 0, my_bodies_count, all_bodies, n);
        for (int i = 0; i < my_bodies_count; i++) {
            my_bodies[i].x += my_bodies[i].vx * dt;
            my_bodies[i].y += my_bodies[i].vy * dt;
            my_bodies[i].z += my_bodies[i].vz * dt;
        }

        // 将每个进程的局部物体数组汇总到全局数组
        MPI_Gather(my_bodies, my_bodies_count, MPI_DOUBLE, all_bodies, my_bodies_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // 计算速度的变化并更新收敛值
        double old_total_velocity_change = total_velocity_change;
        total_velocity_change = compute_total_velocity(my_bodies, 0, n);

        // 输出收敛信息
        if (rank == 0) {
            printf("Convergence at time step %d: %.6f\n", t, total_velocity_change);
        }

        // 如果速度变化小于阈值，则认为系统已经收敛
        if (fabs(old_total_velocity_change - total_velocity_change) < convergence_threshold) {
            if (rank == 0) {
                printf("Converged at time step %d\n", t);
            }
            break;
        }
    }

    end_time = MPI_Wtime();

    // 计算并打印并行版本执行时间，加速比
    if (rank == 0) {
        double parallel_time = end_time - start_time;
        printf("Parallel Execution Time: %.6f seconds\n", parallel_time);
    }

    // 释放内存
    free(all_bodies);
    free(my_bodies);

    MPI_Finalize();

    return 0;
}

