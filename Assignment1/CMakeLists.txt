add_executable(assignment1 qho.c)
target_compile_options(assignment1 PUBLIC -Wall -Wextra -Wpedantic -Werror)
target_link_libraries(assignment1 m)
