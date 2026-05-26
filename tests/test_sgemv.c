#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#define EPSILON 1e-6f
#define NUM_THREADS 4

static int float_equal(float a, float b) {
    return fabsf(a - b) < EPSILON;
}

int test_cblas_sgemv_basic() {
    const int M = 3, N = 2;
    const float alpha = 1.0f, beta = 0.0f;
    float A[] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    float x[] = {1.0f, 2.0f};
    float y[] = {0.0f, 0.0f, 0.0f};
    float expected[] = {5.0f, 11.0f, 17.0f};

    cblas_sgemv(CblasRowMajor, CblasNoTrans, M, N, alpha, A, N, x, 1, beta, y, 1);

    for (int i = 0; i < M; i++) {
        if (!float_equal(y[i], expected[i])) {
            printf("  FAIL: idx=%d got=%f exp=%f\n", i, y[i], expected[i]);
            return 0;
        }
    }
    printf("  PASS: test_cblas_sgemv_basic\n");
    return 1;
}

int test_cblas_sgemv_multithreaded() {
    openblas_set_num_threads(NUM_THREADS);
    const int M = 100, N = 100;
    float *A = malloc(M * N * sizeof(float));
    float *x = malloc(N * sizeof(float));
    float *y = malloc(M * sizeof(float));

    for (int i = 0; i < M * N; i++) A[i] = 1.0f;
    for (int i = 0; i < N; i++) x[i] = 1.0f;
    for (int i = 0; i < M; i++) y[i] = 0.0f;

    cblas_sgemv(CblasRowMajor, CblasNoTrans, M, N, 1.0f, A, N, x, 1, 0.0f, y, 1);

    int passed = 1;
    for (int i = 0; i < M; i++) {
        if (!float_equal(y[i], (float)N)) { passed = 0; break; }
    }
    free(A); free(x); free(y);

    if (passed) {
        printf("  PASS: test_cblas_sgemv_multithreaded (threads=%d)\n", NUM_THREADS);
        return 1;
    }
    printf("  FAIL: test_cblas_sgemv_multithreaded\n");
    return 0;
}