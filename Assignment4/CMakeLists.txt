find_package(MPI)
if (MPI_FOUND)
    set(CMAKE_C_STANDARD 99)
    include_directories(${MPI_INCLUDE_PATH})

    add_executable(assignment4 mt19937ar.c grid.c comms.c makePNG.c rfim.c)
    target_compile_options(assignment4 PUBLIC -Wall -Wextra -Wpedantic ${MPI_COMPILE_FLAGS})
    target_link_options(assignment4 PUBLIC ${MPI_LINK_FLAGS})

    if (NOT CMAKE_BUILD_TYPE MATCHES Debug)
        target_compile_options(assignment4 PUBLIC -Werror)
    endif ()
    
    
    target_link_libraries(assignment4 png m ${MPI_LIBRARIES})
endif ()
