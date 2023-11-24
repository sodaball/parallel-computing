#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 1000
#define tmax 100
#define dt 0.01
#define m 1.0
#define G 6.67430e-11  // 重力常数

// Vector3D是一个包含三个double值的结构体，例如位置或速度
typedef struct {
    double x, y, z;
} Vector3D;

double calculate_distance(Vector3D a, Vector3D b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double dz = b.z - a.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

Vector3D calculate_gravitational_force(Vector3D pos1, Vector3D pos2, double mass1, double mass2) {
    Vector3D force;
    double distance = calculate_distance(pos1, pos2);
    double magnitude = (G * mass1 * mass2) / (distance * distance);

    force.x = magnitude * (pos2.x - pos1.x) / distance;
    force.y = magnitude * (pos2.y - pos1.y) / distance;
    force.z = magnitude * (pos2.z - pos1.z) / distance;

    return force;
}

void update_positions_velocities(Vector3D *x, Vector3D *v, Vector3D *xnew, Vector3D *vnew, int start, int end) {
    for (int i = start; i < end; i++) {
        vnew[i] = v[i];  // 用当前速度初始化vnew

        for (int j = 0; j < n; j++) {
            if (i != j) {
                // 计算本进程负责的质点收到的力，force为结构体，包含三个double值
                Vector3D force = calculate_gravitational_force(x[i], x[j], m, m);
                // 计算力之后更新速度
                vnew[i].x += force.x * dt / m;
                vnew[i].y += force.y * dt / m;
                vnew[i].z += force.z * dt / m;
            }
        }
        // 计算速度之后更新位置
        xnew[i].x = x[i].x + vnew[i].x * dt;
        xnew[i].y = x[i].y + vnew[i].y * dt;
        xnew[i].z = x[i].z + vnew[i].z * dt;
    }
}

void random_initialize(Vector3D *array, int size) {
    for (int i = 0; i < size; i++) {
        array[i].x = (double)rand() / RAND_MAX;  // 0到1之间的随机值
        array[i].y = (double)rand() / RAND_MAX;
        array[i].z = (double)rand() / RAND_MAX;
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    double start, end;

    start = MPI_Wtime();

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_n = n / size;  // 假设n可以被进程数整除

    // 为位置和速度分配内存
    Vector3D *x = (Vector3D *)malloc(n * sizeof(Vector3D)); // 每个质点的位置
    Vector3D *v = (Vector3D *)malloc(n * sizeof(Vector3D)); // 每个质点的速度
    Vector3D *xnew = (Vector3D *)malloc(n * sizeof(Vector3D));  // 每个质点位置的更新
    Vector3D *vnew = (Vector3D *)malloc(n * sizeof(Vector3D));  // 每个质点速度的更新

    // 初始化位置和速度
    random_initialize(x, n);
    random_initialize(v, n);

    // 主时间步进循环
    for (int t = 0; t < tmax; t++) {
        update_positions_velocities(x, v, xnew, vnew, rank * local_n, (rank + 1) * local_n);

        // 同步所有进程的数据，包括位置和速度
        MPI_Allgather(xnew + rank * local_n, local_n, MPI_DOUBLE, x, local_n, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(vnew + rank * local_n, local_n, MPI_DOUBLE, v, local_n, MPI_DOUBLE, MPI_COMM_WORLD);

        // 输出信息（例如，每10步一次）
        if (t % 100 == 0) {
            for (int i = 0; i < n; i++) {
                printf("Process %d, Particle %d: Position (%f, %f, %f), Velocity (%f, %f, %f)\n",
                   rank, i, x[i].x, x[i].y, x[i].z, v[i].x, v[i].y, v[i].z);
            }
        }
    }

    // 清理内存
    free(x);
    free(v);
    free(xnew);
    free(vnew);

    end = MPI_Wtime();

    if (rank == 0) {
        printf("进程数: %d\n", size);
        printf("质点数: %d\n", n);
        printf("时间步数: %d\n", tmax);
        printf("总时间: %f 秒.\n", end - start);
    }

    MPI_Finalize();
    return 0;
}
