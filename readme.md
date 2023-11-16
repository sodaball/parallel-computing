### 软硬件环境：

![image-20231114111421488](C:\Users\Dachau\AppData\Roaming\Typora\typora-user-images\image-20231114111421488.png)

内核：Linux

网络节点上的主机名：kunpeng66.localdomin

内核发行号：4.19.90

内核版本：#4 SMP Tue Feb 15 22:54:24 CST 2022

硬件架构：aarch64

操作系统：GNU/Linux



### 题目一

#### pthread

编译：

```bash
gcc -o riemann_zeta_pthread riemann_zeta_pthread.c -lpthread -lm
```

运行：

```bash
./riemann_zeta_pthread num_threads num_s num_k
```

结果：

![image-20231114141524614](C:\Users\Dachau\AppData\Roaming\Typora\typora-user-images\image-20231114141524614.png)

上图是x和k不变的情况下，只改变线程数量





### 题目二

#### OpenMP

编译：

```
gcc -fopenmp -o ./nbody_openmp ./nbody_openmp.c  -lm
```



运行：

```
./nbody_openmp
```



问题1：

在输出信息中，出现了"nan"，这表示在计算中出现了非数值（Not-a-Number）。这通常是由于物体之间的距离过小而导致的数值不稳定性，导致计算引力的分母变为零，进而导致除以零的情况。为了避免这种情况，可以在计算前添加一些检查和处理。



问题2：

加速比有时候大于1有时候小于1

串行和并行执行两次的初始值不一样，应该设置为一样的



结果：（**输出时间步为1000，可以调整**）

![image-20231103173551059](C:\Users\Dachau\AppData\Roaming\Typora\typora-user-images\image-20231103173551059.png)



### 题目三 

#### OpenMP

编译：

```bash
gcc -fopenmp -o ./kmeans_iris_openmp ./kmeans_iris_openmp.c -lm
```

`#pragma omp parallel for`是OpenMP中的一个指令，表示接下来的for循环将被多线程执行，另外每次循环之间不能有关系

运行：

```bash
./kmeans_iris_openmp
```



问题1：

/usr/bin/ld: /tmp/cc9Mgky3.o: in function euclidean_distance:
kmeans.c:(.text+0x80): undefined reference to sqrt

解决：

编译时加 `-lm`



问题2：

多线程比单线程还慢，可能是线程管理开销比较大，但是尝试增加NUM_POINTS或MAX_ITERATIONS的值

![image-20231114112203951](C:\Users\Dachau\AppData\Roaming\Typora\typora-user-images\image-20231114112203951.png)
