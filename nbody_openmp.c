#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define G 6.67430e-11  // 万有引力常数
#define tmax 1000       // 最大时间步数
#define dt 1.0         // 时间步长
#define convergence_threshold 1e-6  // 收敛阈值

typedef struct {
    double mass;
    double x, y, z;
    double vx, vy, vz;
} Body;

void initialize_bodies(Body* bodies, int num_bodies) {
    // 设置相同的随机种子
    srand(1234);  // 可以选择任何非负整数

    #pragma omp parallel for
    for (int i = 0; i < num_bodies; i++) {
        bodies[i].mass = 1e24;  // 随机设置质量
        bodies[i].x = rand() / (double)RAND_MAX;  // 随机初始位置
        bodies[i].y = rand() / (double)RAND_MAX;
        bodies[i].z = rand() / (double)RAND_MAX;
        bodies[i].vx = 0.0;
        bodies[i].vy = 0.0;
        bodies[i].vz = 0.0;
    }
}

void compute_force(Body* my_bodies, int my_start, int my_end, Body* all_bodies, int n) {
    double softening = 1e-6;  // 软化参数

    #pragma omp parallel for
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
    // 设置线程数
    omp_set_num_threads(atoi(argv[1]));

    int n = 1000; // N体问题中的物体数量

    // 创建全局物体数组
    Body* all_bodies = (Body*)malloc(sizeof(Body) * n);
    // 初始化所有物体的质量、位置和速度信息
    initialize_bodies(all_bodies, n);

    // 创建每个进程负责的局部物体数组
    Body* my_bodies = (Body*)malloc(sizeof(Body) * n);
    // 复制全局数组到每个进程的局部数组
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        my_bodies[i] = all_bodies[i];
    }

    double total_velocity_change = 1.0;  // 初始值，保证至少一次循环
    double start_time, end_time;

    start_time = omp_get_wtime();

    for (int t = 0; t < tmax && total_velocity_change > convergence_threshold; t++) {
        // 计算受力和更新位置
        compute_force(my_bodies, 0, n, all_bodies, n);
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            my_bodies[i].x += my_bodies[i].vx * dt;
            my_bodies[i].y += my_bodies[i].vy * dt;
            my_bodies[i].z += my_bodies[i].vz * dt;
        }

        // 计算速度的变化并更新收敛值
        double old_total_velocity_change = total_velocity_change;
        total_velocity_change = compute_total_velocity(my_bodies, 0, n);

        // 输出每个时间步结束时的部分位置信息
        // if (t % 1000 == 0) {
        //     #pragma omp parallel for
        //     for (int i = 0; i < n; i++) {
        //         printf("Body %d: x=%.2f, y=%.2f, z=%.2f\n", i, my_bodies[i].x, my_bodies[i].y, my_bodies[i].z);
        //     }
        // }

        // 输出收敛信息
        printf("Convergence at time step %d: %.6f\n", t, total_velocity_change);

        // 如果速度变化小于阈值，则认为系统已经收敛
        if (fabs(old_total_velocity_change - total_velocity_change) < convergence_threshold) {
            printf("Converged at time step %d\n", t);
            break;
        }
    }

    end_time = omp_get_wtime();

    // 计算并打印串行版本执行时间
    double serial_time = end_time - start_time;

    start_time = omp_get_wtime();

    // 执行并行版本
    // initialize_bodies(all_bodies, n);
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        my_bodies[i] = all_bodies[i];
    }

    total_velocity_change = 1.0;  // 重置收敛值
    for (int t = 0; t < tmax && total_velocity_change > convergence_threshold; t++) {
        // 计算受力和更新位置
        compute_force(my_bodies, 0, n, all_bodies, n);
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            my_bodies[i].x += my_bodies[i].vx * dt;
            my_bodies[i].y += my_bodies[i].vy * dt;
            my_bodies[i].z += my_bodies[i].vz * dt;
        }

        // 计算速度的变化并更新收敛值
        double old_total_velocity_change = total_velocity_change;
        total_velocity_change = compute_total_velocity(my_bodies, 0, n);

        // // 输出每个时间步结束时的部分位置信息
        // if (t % 1000 == 0) {
        //     #pragma omp parallel for
        //     for (int i = 0; i < n; i++) {
        //         printf("Body %d: x=%.2f, y=%.2f, z=%.2f\n", i, my_bodies[i].x, my_bodies[i].y, my_bodies[i].z);
        //     }
        // }

        // 输出收敛信息
        printf("Convergence at time step %d: %.6f\n", t, total_velocity_change);

        // 如果速度变化小于阈值，则认为系统已经收敛
        if (fabs(old_total_velocity_change - total_velocity_change) < convergence_threshold) {
            printf("Converged at time step %d\n", t);
            break;
        }
    }

    end_time = omp_get_wtime();

    // 计算并打印并行版本执行时间，加速比
    printf("Serial Execution Time: %.6f seconds\n", serial_time);
    double parallel_time = end_time - start_time;
    printf("Parallel Execution Time: %.6f seconds\n", parallel_time);

    double speedup = serial_time / parallel_time;
    printf("Speedup: %.2f\n", speedup);

    // 释放内存
    free(all_bodies);
    free(my_bodies);

    // 输出进程数、物体个数和迭代次数
    printf("Number of processes: %d\n", omp_get_max_threads());
    printf("Number of bodies: %d\n", n);
    printf("Number of iterations: %d\n", tmax);

    return 0;
}
