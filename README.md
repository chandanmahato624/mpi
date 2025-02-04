

1. WAPP in C/C++ to compute PI by the Monte Carlo method, testing whether points in the unit square are in the unit
circle

```javascript
#include <mpi.h>
#include <stdio.h>
#incl#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
void bubble_sort(int arr[], int n) {
int i, j;
for (i = 0; i < n-1; i++) {
for (j = 0; j < n-i-1; j++) {
if (arr[j] > arr[j+1]) {
int temp = arr[j];
arr[j] = arr[j+1];
arr[j+1] = temp; }}}}
int main(int argc, char *argv[]) {
int pid, np;
int *arr = NULL;
int n = 10;
int i;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &pid);
MPI_Comm_size(MPI_COMM_WORLD, &np);
if (pid == 0) {
arr = (int *)malloc(n * sizeof(int));
for (i = 0; i < n; i++) {
arr[i] = rand() % 100; }
printf("Original array:\n");
for (i = 0; i < n; i++) {
printf("%d ", arr[i]);}
printf("\n"); }
MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
int local_n = n / np;
int *local_arr = (int *)malloc(local_n * sizeof(int));
MPI_Scatter(arr, local_n, MPI_INT, local_arr, local_n, MPI_INT, 0, 
MPI_COMM_WORLD);
bubble_sort(local_arr, local_n);
MPI_Gather(local_arr, local_n, MPI_INT, arr, local_n, MPI_INT, 0, 
MPI_COMM_WORLD);
if (pid == 0) {
bubble_sort(arr, n);
printf("Sorted array:\n");
for (i = 0; i < n; i++) {
printf("%d ", arr[i]);}
printf("\n");}
MPI_Finalize();
return 0;}
```

2. WAPP in C/C++ to compute the sum of elements in an array using scatter and gather operation using MPI
for parallel execution.

```javascript
https://colab.research.google.com/drive/1Bjb5zbwpRR7pEJT_NE70Q5sC2hAcvBNg?usp=sharing#scrollTo=E22r8nUtyExO
```


3. WAPP in C/C++ to implement parallel linear search using MPI for parallel execution.

```javascript
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 16
int main(int argc, char* argv[]) {
int rank, size, i;
int a[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
int x = 10; // Element to search for
int local_index = -1, global_index = -1;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
int elements_per_process = N / size;
int sub_a[elements_per_process];
MPI_Scatter(a, elements_per_process, MPI_INT, sub_a, elements_per_process, MPI_INT, 0,
MPI_COMM_WORLD);
for (i = 0; i < elements_per_process; i++) {
if (sub_a[i] == x) {
local_index = rank * elements_per_process + i;
break;
}
}
MPI_Allreduce(&local_index, &global_index, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
if (rank == 0) {
if (global_index != -1) {
printf("Element found at index %d\n", global_index);
} else {
printf("Element not found\n");
}
}
MPI_Finalize();
return 0;
}
```



4. WAPP in C/C++ to find the largest value among the elements stored in an array using scatter and gather operation
using MPI for parallel execution.

```javascript
#include <mpi.h>
#include <stdio.h>
#define N 16
int main(int argc, char* argv[]) {
int rank, size, i;
int a[N] = {1, 2, 3, 4, 15, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
int local_max, global_max;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
int elements_per_process = N / size;
int sub_a[elements_per_process];
MPI_Scatter(a, elements_per_process, MPI_INT, sub_a, elements_per_process, MPI_INT, 0,
MPI_COMM_WORLD);
local_max = sub_a[0];
for (i = 1; i < elements_per_process; i++) {
if (sub_a[i] > local_max) {
local_max = sub_a[i];
}
}
MPI_Reduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
if (rank == 0) {
printf("Maximum value = %d\n", global_max);
}
MPI_Finalize();
return
```


5. .WAPP in C/C++ using MPI to compute the sum of all odd numbers in an array using broadcast.

```javascript
#include <mpi.h>
#include <stdio.h>
#define N 16
int main(int argc, char* argv[]) {
int rank, size, i;
int a[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
int sub_sum = 0, total_sum = 0;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
int elements_per_process = N / size;
int sub_a[elements_per_process];
MPI_Scatter(a, elements_per_process, MPI_INT, sub_a, elements_per_process, MPI_INT, 0,
MPI_COMM_WORLD);
for (i = 0; i < elements_per_process; i++) {
if (sub_a[i] % 2 != 0) {
sub_sum += sub_a[i];
}
}
MPI_Reduce(&sub_sum, &total_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
if (rank == 0) {
printf("Total sum of odd numbers = %d\n", total_sum);
}
MPI_Finalize();
return 0;
}
```

6. Mat Add

```javascript
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 4
void print_matrix(int mat[N][N]) {
for (int i = 0; i < N; i++) {
for (int j = 0; j < N; j++) {
printf("%d ", mat[i][j]);
}
printf("\n");
}
}
int main(int argc, char* argv[]) {
int rank, size, i, j, k;
int a[N][N] = {
{1, 2, 3, 4},
{5, 6, 7, 8},
{9, 10, 11, 12},
{13, 14, 15, 16}
};
int b[N][N] = {
{1, 0, 0, 0},
{0, 1, 0, 0},
{0, 0, 1, 0},
{0, 0, 0, 1}
};
int c[N][N] = {0};
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
if (N % size != 0) {
if (rank == 0) {
printf("Matrix size must be divisible by the number of processes.\n");
}
MPI_Finalize();
return -1;
}
int rows_per_process = N / size;
int sub_a[rows_per_process][N], sub_c[rows_per_process][N];
MPI_Scatter(a, rows_per_process * N, MPI_INT, sub_a, rows_per_process * N, MPI_INT, 0,
MPI_COMM_WORLD);
MPI_Bcast(b, N * N, MPI_INT, 0, MPI_COMM_WORLD);
for (i = 0; i < rows_per_process; i++) {
for (j = 0; j < N; j++) {
sub_c[i][j] = 0;
for (k = 0; k < N; k++) {
sub_c[i][j] += sub_a[i][k] * b[k][j];
}
}
}
MPI_Gather(sub_c, rows_per_process * N, MPI_INT, c, rows_per_process * N, MPI_INT, 0,
MPI_COMM_WORLD);
if (rank == 0) {
printf("Result matrix:\n");
print_matrix(c);
}
MPI_Finalize();
return 0;
}
```


7. Mat sum

```javascript
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 4
void print_matrix(int mat[N][N]) {
for (int i = 0; i < N; i++) {
for (int j = 0; j < N; j++) {
printf("%d ", mat[i][j]);
}
printf("\n");
}
}
int main(int argc, char* argv[]) {
int rank, size, i, j;
int a[N][N] = {
{1, 2, 3, 4},
{5, 6, 7, 8},
{9, 10, 11, 12},
{13, 14, 15, 16}
};
int b[N][N] = {
{16, 15, 14, 13},
{12, 11, 10, 9},
{8, 7, 6, 5},
{4, 3, 2, 1}
};
int c[N][N] = {0};
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
if (N % size != 0) {
if (rank == 0) {
printf("Matrix size must be divisible by the number of processes.\n");
}
MPI_Finalize();
return -1;
}
int rows_per_process = N / size;
int sub_a[rows_per_process][N], sub_b[rows_per_process][N], sub_c[rows_per_process][N];
MPI_Scatter(a, rows_per_process * N, MPI_INT, sub_a, rows_per_process * N, MPI_INT, 0,
MPI_COMM_WORLD);
MPI_Scatter(b, rows_per_process * N, MPI_INT, sub_b, rows_per_process * N, MPI_INT, 0,
MPI_COMM_WORLD);
for (i = 0; i < rows_per_process; i++) {
for (j = 0; j < N; j++) {
sub_c[i][j] = sub_a[i][j] + sub_b[i][j];
}
}
MPI_Gather(sub_c, rows_per_process * N, MPI_INT, c, rows_per_process * N, MPI_INT, 0,
MPI_COMM_WORLD);
if (rank == 0) {
printf("Result matrix:\n");
print_matrix(c);
}
MPI_Finalize();
return 0;
}
```


####1. WAPP in C/C++ to compute PI by the Monte Carlo method, testing whether points in the unit square are in the unit
circle

```javascript
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4

void print_matrix(int mat[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", mat[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]) {
    int rank, size, i, j, k;
    int a[N][N] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12},
        {13, 14, 15, 16}
    };
    int b[N][N] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    int c[N][N] = {0};

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (N % size != 0) {
        if (rank == 0) {
            printf("Matrix size must be divisible by the number of processes.\n");
        }
        MPI_Finalize();
        return -1;
    }

    int rows_per_process = N / size;
    int sub_a[rows_per_process][N], sub_c[rows_per_process][N];

    MPI_Scatter(a, rows_per_process * N, MPI_INT, sub_a, rows_per_process * N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    for (i = 0; i < rows_per_process; i++) {
        for (j = 0; j < N; j++) {
            sub_c[i][j] = 0;
            for (k = 0; k < N; k++) {
                sub_c[i][j] += sub_a[i][k] * b[k][j];
            }
        }
    }

    MPI_Gather(sub_c, rows_per_process * N, MPI_INT, c, rows_per_process * N, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Result matrix:\n");
        print_matrix(c);
    }

    MPI_Finalize();
    return 0;
}

```
