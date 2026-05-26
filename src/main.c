#include <stdio.h>

// Объявления тестов
int test_cblas_sgemv_basic();
int test_cblas_sgemv_multithreaded();
int test_cblas_dgemv_basic();
int test_cblas_sgemm_basic();
int test_cblas_dgemm_basic();

int main() {
    int passed = 0, total = 0;
    
    printf("=== CBLAS Interface Tests ===\n\n");
    
    // SGEMV tests
    total++; passed += test_cblas_sgemv_basic();
    total++; passed += test_cblas_sgemv_multithreaded();
    
    // DGEMV tests
    total++; passed += test_cblas_dgemv_basic();
    
    // SGEMM tests
    total++; passed += test_cblas_sgemm_basic();
    
    // DGEMM tests
    total++; passed += test_cblas_dgemm_basic();
    
    printf("\n=== Results: %d/%d tests passed ===\n", passed, total);
    return (passed == total) ? 0 : 1;
}