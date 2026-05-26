#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#define EPSILON 1e-12

static int double_equal(double a, double b) {
    return fabs(a - b) < EPSILON;
}

int test_cblas_dgemv_basic() {
    const int M = 2, N = 3;
    const double alpha = 2.0, beta = 0.5;
    double A[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double x[] = {1.0, 2.0, 3.0};
    double y[] = {1.0, 1.0};
    double expected[] = {28.5, 64.5};

    cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, alpha, A, N, x, 1, beta, y, 1);

    for (int i = 0; i < M; i++) {
        if (!double_equal(y[i], expected[i])) {
            printf("  FAIL: idx=%d got=%f exp=%f\n", i, y[i], expected[i]);
            return 0;
        }
    }
    printf("  PASS: test_cblas_dgemv_basic\n");
    return 1;
}