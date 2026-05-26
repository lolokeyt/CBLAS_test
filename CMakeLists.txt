cmake_minimum_required(VERSION 3.10)
project(CBLAS_Test C)

set(CMAKE_C_STANDARD 11)

# Поиск OpenBLAS
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(OPENBLAS openblas)
endif()

if(OPENBLAS_FOUND)
    include_directories(${OPENBLAS_INCLUDE_DIRS})
    link_directories(${OPENBLAS_LIBRARY_DIRS})
else()
    # fallback для Windows / ручная установка
    include_directories("C:/OpenBLAS/include")
    link_directories("C:/OpenBLAS/lib")
endif()

# Исходные файлы тестов
set(TEST_SOURCES
    src/main.c
    tests/test_sgemv.c
    tests/test_dgemv.c
    tests/test_sgemm.c
    tests/test_dgemm.c
)

add_executable(test_cblas ${TEST_SOURCES})

# Линковка с OpenBLAS
if(WIN32)
    target_link_libraries(test_cblas openblas)
else()
    target_link_libraries(test_cblas openblas pthread m)
endif()

# Тесты
enable_testing()
add_test(NAME cblas_tests COMMAND test_cblas)