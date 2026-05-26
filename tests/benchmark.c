#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "cblas.h"

// Объявления ваших функций
void my_ssymm(const enum CBLAS_ORDER, const enum CBLAS_SIDE,
              const enum CBLAS_UPLO, const int, const int,
              const float, const float*, const int,
              const float*, const int, const float,
              float*, const int);

void my_dsymm(const enum CBLAS_ORDER, const enum CBLAS_SIDE,
              const enum CBLAS_UPLO, const int, const int,
              const double, const double*, const int,
              const double*, const int, const double,
              double*, const int);

void init_matrix_float(float *A, int rows, int cols) {
    for (int i = 0; i < rows * cols; i++) {
        A[i] = (float)rand() / RAND_MAX * 2.0 - 1.0;
    }
}

void init_matrix_double(double *A, int rows, int cols) {
    for (int i = 0; i < rows * cols; i++) {
        A[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }
}

double get_time() {
    return omp_get_wtime();
}

void benchmark_ssymm(int size, int iterations, int threads) {
    omp_set_num_threads(threads);
    
    int M = size, N = size;
    float *A = malloc(M * M * sizeof(float));
    float *B = malloc(M * N * sizeof(float));
    float *C_my = malloc(M * N * sizeof(float));
    float *C_blas = malloc(M * N * sizeof(float));
    
    init_matrix_float(A, M, M);
    init_matrix_float(B, M, N);
    
    printf("\nSSYMM\nРазмер массива: %dx%d\n", size, size);
    printf("Количество итераций: %d\n", iterations);
    printf("Потоки: %d\n\n", threads);
    printf("| Запуск | Реализация (с) | OpenBLAS (с) | Производительность %% | Отношение |\n");
    printf("|--------|----------------|--------------|---------------------|-----------|\n");
    
    double total_ratio = 0.0;
    
    for (int iter = 0; iter < iterations; iter++) {
        // Ваша реализация
        double start_my = get_time();
        my_ssymm(CblasRowMajor, CblasLeft, CblasUpper, M, N,
                 1.0f, A, M, B, N, 0.0f, C_my, N);
        double time_my = get_time() - start_my;
        
        // OpenBLAS
        double start_blas = get_time();
        cblas_ssymm(CblasRowMajor, CblasLeft, CblasUpper, M, N,
                    1.0f, A, M, B, N, 0.0f, C_blas, N);
        double time_blas = get_time() - start_blas;
        
        double ratio = time_my / time_blas;
        double perf_percent = (time_blas / time_my) * 100.0;
        total_ratio += ratio;
        
        printf("| %6d | %14.6f | %12.6f | %19.2f | %9.2f |\n",
               iter + 1, time_my, time_blas, perf_percent, ratio);
    }
    
    double geom_mean = pow(total_ratio / iterations, 1.0);
    printf("\nСреднее геометрическое производительности: %.2f%%\n", 100.0 / geom_mean);
    printf("Ускорение относительно OpenBLAS: %.2fx\n", 1.0 / geom_mean);
    
    free(A); free(B); free(C_my); free(C_blas);
}

void benchmark_dsymm(int size, int iterations, int threads) {
    omp_set_num_threads(threads);
    
    int M = size, N = size;
    double *A = malloc(M * M * sizeof(double));
    double *B = malloc(M * N * sizeof(double));
    double *C_my = malloc(M * N * sizeof(double));
    double *C_blas = malloc(M * N * sizeof(double));
    
    init_matrix_double(A, M, M);
    init_matrix_double(B, M, N);
    
    printf("\nDSYMM\nРазмер массива: %dx%d\n", size, size);
    printf("Количество итераций: %d\n", iterations);
    printf("Потоки: %d\n\n", threads);
    printf("| Запуск | Реализация (с) | OpenBLAS (с) | Производительность %% | Отношение |\n");
    printf("|--------|----------------|--------------|---------------------|-----------|\n");
    
    double total_ratio = 0.0;
    
    for (int iter = 0; iter < iterations; iter++) {
        double start_my = get_time();
        my_dsymm(CblasRowMajor, CblasLeft, CblasUpper, M, N,
                 1.0, A, M, B, N, 0.0, C_my, N);
        double time_my = get_time() - start_my;
        
        double start_blas = get_time();
        cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, M, N,
                    1.0, A, M, B, N, 0.0, C_blas, N);
        double time_blas = get_time() - start_blas;
        
        double ratio = time_my / time_blas;
        double perf_percent = (time_blas / time_my) * 100.0;
        total_ratio += ratio;
        
        printf("| %6d | %14.6f | %12.6f | %19.2f | %9.2f |\n",
               iter + 1, time_my, time_blas, perf_percent, ratio);
    }
    
    double geom_mean = pow(total_ratio / iterations, 1.0);
    printf("\nСреднее геометрическое производительности: %.2f%%\n", 100.0 / geom_mean);
    printf("Ускорение относительно OpenBLAS: %.2fx\n", 1.0 / geom_mean);
    
    free(A); free(B); free(C_my); free(C_blas);
}

int main() {
    int size = 2500;  // Подберите размер для 1+ минуты выполнения
    int iterations = 10;
    
    printf("=== Бенчмарк SYMM ===\n");
    
    int threads[] = {1, 2, 4, 8, 16};
    int num_configs = sizeof(threads) / sizeof(threads[0]);
    
    for (int i = 0; i < num_configs; i++) {
        benchmark_ssymm(size, iterations, threads[i]);
        benchmark_dsymm(size, iterations, threads[i]);
    }
    
    return 0;
}