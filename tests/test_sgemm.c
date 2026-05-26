#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#define EPSILON 1e-6f

static int float_equal(float a, float b) {
    return fabsf(a - b) < EPSILON;
}

int test_cblas_sgemm_basic() {
    const int M = 2, N = 3, K = 4;
    const float alpha = 1.0f, beta = 0.0f;
    float A[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    float B[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    float C[6] = {0};
    float expected[6] = {70, 80, 90, 158, 184, 210};

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                M, N, K, alpha, A, K, B, N, beta, C, N);

    for (int i = 0; i < M * N; i++) {
        if (!float_equal(C[i], expected[i])) {
            printf("  FAIL: idx=%d got=%f exp=%f\n", i, C[i], expected[i]);
            return 0;
        }
    }
    printf("  PASS: test_cblas_sgemm_basic\n");
    return 1;
}