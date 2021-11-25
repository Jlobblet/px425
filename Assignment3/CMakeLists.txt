set(CMAKE_C_STANDARD 99)

find_package(OpenMP)
list(APPEND CompileOptions -Wall -Wextra -Wpedantic -Wno-unknown-pragmas -Werror)
list(APPEND Libraries pthread gfortran rt)

add_executable(assignment3serial percolation.c mt19937ar.c cluster_search.f90)
target_compile_options(assignment3serial PUBLIC ${CompileOptions})
target_link_libraries(assignment3serial ${Libraries})

if (OpenMP_C_FOUND)
    add_executable(assignment3omp percolation.c mt19937ar.c cluster_search.f90)
    target_compile_options(assignment3omp PUBLIC ${CompileOptions} ${OpenMP_C_FLAGS})
    target_link_libraries(assignment3omp ${Libraries} OpenMP::OpenMP_C)

    add_executable(assignment3.1omp percolation1.c mt19937ar.c cluster_search.f90)
    target_compile_options(assignment3.1omp PUBLIC ${CompileOptions} ${OpenMP_C_FLAGS})
    target_link_libraries(assignment3.1omp ${Libraries} OpenMP::OpenMP_C)

    add_executable(assignment3.2omp percolation2.c mt19937ar.c cluster_search.f90)
    target_compile_options(assignment3.2omp PUBLIC ${CompileOptions} ${OpenMP_C_FLAGS})
    target_link_libraries(assignment3.2omp ${Libraries} OpenMP::OpenMP_C)
endif()
