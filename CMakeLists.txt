cmake_minimum_required(VERSION 3.10)
project(CBLAS_Test C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

# Пути к OpenBLAS
include_directories(/opt/OpenBLAS/include)
link_directories(/opt/OpenBLAS/lib)

# Ваша реализация
add_library(my_symm src/my_symm.c)

# Бенчмарк
add_executable(benchmark tests/benchmark.c)
target_link_libraries(benchmark my_symm openblas m)