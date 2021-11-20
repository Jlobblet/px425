set(CMAKE_C_STANDARD 99)

list(APPEND CompileOptions -Wall -Wextra -Wpedantic -Wno-unknown-pragmas -Werror)
list(APPEND Libraries pthread gfortran rt)

add_executable(assignment3serial percolation.c mt19937ar.c cluster_search.f90)
target_compile_options(assignment3serial PUBLIC ${CompileOptions})
target_link_libraries(assignment3serial ${Libraries})

add_executable(assignment3omp percolation.c mt19937ar.c cluster_search.f90)
target_compile_options(assignment3omp PUBLIC ${CompileOptions} ${OpenMP_C_FLAGS})
target_link_libraries(assignment3omp ${Libraries})
target_link_options(assignment3omp PUBLIC ${OpenMP_EXE_LINKER_FLAGS})
