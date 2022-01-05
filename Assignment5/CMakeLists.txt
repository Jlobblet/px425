find_package(MPI)
find_package(OpenMP)

if (MPI_FOUND AND OpenMP_C_FOUND)
    set(CMAKE_C_STANDARD 99)
    include_directories(${MPI_INCLUDE_PATH})

    # Sanitisers
    list(APPEND AsanOptions -fno-omit-frame-pointer -fsanitize=address)
    list(APPEND LsanOptions -fsanitize=leak)
    list(APPEND TsanOptions -fsanitize=thread)
    list(APPEND UbsanOptions -fsanitize=undefined)
    list(APPEND MsanOptions -fsanitize=memory)

    list(APPEND CompileOptions -Wall -Wextra -Wpedantic -Wno-unused-parameter -mavx2 -march=native ${OpenMP_C_FLAGS} ${MPI_COMPILE_FLAGS})

    list(APPEND LinkOptions ${MPI_LINK_FLAGS})
    list(APPEND LinkLibraries ${MPI_LIBRARIES} OpenMP::OpenMP_C m)

    if (CMAKE_BUILD_TYPE MATCHES DEBUG)
        list(APPEND CompileOptions ${TsanOptions})
        list(APPEND LinkOptions ${TsanOptions})
    else ()
        list(APPEND CompileOptions -Werror)
    endif ()

    add_executable(assignment5 wifi.c mt19937ar.c input.c output.c comms.c types.c)
    target_compile_options(assignment5 PUBLIC ${CompileOptions})
    target_link_options(assignment5 PUBLIC ${LinkOptions})
    target_link_libraries(assignment5 ${LinkLibraries})
endif ()
