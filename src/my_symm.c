#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"

// Реализация SSYMM (одинарная точность)
void my_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
              const enum CBLAS_UPLO Uplo, const int M, const int N,
              const float alpha, const float *A, const int lda,
              const float *B, const int ldb, const float beta,
              float *C, const int ldc) {
    if (M == 0 || N == 0) return;
    if (alpha == 0.0f && beta == 1.0f) return;
    
    if (alpha == 0.0f) {
        if (beta == 0.0f) {
            for (int j = 0; j < N; j++)
                for (int i = 0; i < M; i++)
                    C[i + j * ldc] = 0.0f;
        } else {
            for (int j = 0; j < N; j++)
                for (int i = 0; i < M; i++)
                    C[i + j * ldc] = beta * C[i + j * ldc];
        }
        return;
    }
    
    if (Order == CblasRowMajor) {
        if (Side == CblasLeft) {
            if (Uplo == CblasUpper) {
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < M; i++) {
                        float temp = 0.0f;
                        for (int k = 0; k < i; k++)
                            temp += A[i * lda + k] * B[k + j * ldb];
                        for (int k = i; k < M; k++)
                            temp += A[k * lda + i] * B[k + j * ldb];
                        C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
                    }
                }
            } else {
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < M; i++) {
                        float temp = 0.0f;
                        for (int k = 0; k <= i; k++)
                            temp += A[i * lda + k] * B[k + j * ldb];
                        for (int k = i + 1; k < M; k++)
                            temp += A[k * lda + i] * B[k + j * ldb];
                        C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
                    }
                }
            }
        } else {
            if (Uplo == CblasUpper) {
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < M; i++) {
                        float temp = 0.0f;
                        for (int k = 0; k < j; k++)
                            temp += B[i + k * ldb] * A[k * lda + j];
                        for (int k = j; k < N; k++)
                            temp += B[i + k * ldb] * A[j * lda + k];
                        C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
                    }
                }
            } else {
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < M; i++) {
                        float temp = 0.0f;
                        for (int k = 0; k <= j; k++)
                            temp += B[i + k * ldb] * A[k * lda + j];
                        for (int k = j + 1; k < N; k++)
                            temp += B[i + k * ldb] * A[j * lda + k];
                        C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
                    }
                }
            }
        }
    }
}

// Реализация DSYMM (двойная точность)
void my_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
              const enum CBLAS_UPLO Uplo, const int M, const int N,
              const double alpha, const double *A, const int lda,
              const double *B, const int ldb, const double beta,
              double *C, const int ldc) {
    if (M == 0 || N == 0) return;
    if (alpha == 0.0 && beta == 1.0) return;
    
    if (alpha == 0.0) {
        if (beta == 0.0) {
            for (int j = 0; j < N; j++)
                for (int i = 0; i < M; i++)
                    C[i + j * ldc] = 0.0;
        } else {
            for (int j = 0; j < N; j++)
                for (int i = 0; i < M; i++)
                    C[i + j * ldc] = beta * C[i + j * ldc];
        }
        return;
    }
    
    if (Order == CblasRowMajor) {
        if (Side == CblasLeft) {
            if (Uplo == CblasUpper) {
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < M; i++) {
                        double temp = 0.0;
                        for (int k = 0; k < i; k++)
                            temp += A[i * lda + k] * B[k + j * ldb];
                        for (int k = i; k < M; k++)
                            temp += A[k * lda + i] * B[k + j * ldb];
                        C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
                    }
                }
            } else {
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < M; i++) {
                        double temp = 0.0;
                        for (int k = 0; k <= i; k++)
                            temp += A[i * lda + k] * B[k + j * ldb];
                        for (int k = i + 1; k < M; k++)
                            temp += A[k * lda + i] * B[k + j * ldb];
                        C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
                    }
                }
            }
        }
    }
}

