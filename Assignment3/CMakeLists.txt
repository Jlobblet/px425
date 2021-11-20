set(CMAKE_C_STANDARD 99)

add_executable(assignment3 percolation.c mt19937ar.c cluster_search.f90)
target_compile_options(assignment3 PUBLIC -Wall -Wextra -Wpedantic)
target_link_libraries(assignment3 pthread gfortran rt)
