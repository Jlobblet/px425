add_compile_options(-Wall -Wextra -Wpedantic -Werror)

add_executable(assignment1 qho.c)
target_link_libraries(assignment1 m)
