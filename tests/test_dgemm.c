#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#define EPSILON 1e-12

static int double_equal(double a, double b) {
    return fabs(a - b) < EPSILON;
}

int test_cblas_dgemm_basic() {
    const int M = 2, N = 2, K = 2;
    const double alpha = 1.0, beta = 0.0;
    double A[4] = {1, 2, 3, 4};
    double B[4] = {5, 6, 7, 8};
    double C[4] = {0};
    double expected[4] = {19, 22, 43, 50};

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, alpha, A, K, B, N, beta, C, N);

    for (int i = 0; i < M * N; i++) {
        if (!double_equal(C[i], expected[i])) {
            printf("  FAIL: idx=%d\n", i);
            return 0;
        }
    }
    printf("  PASS: test_cblas_dgemm_basic\n");
    return 1;
}