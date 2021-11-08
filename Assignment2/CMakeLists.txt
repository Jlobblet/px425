set(CMAKE_C_STANDARD 99)
add_compile_options(-Wall -Wextra -Wpedantic)

add_executable(assignment2 burgers.c makePNG.c mt19937ar.c)
target_link_libraries(assignment2 m png)
add_custom_target(assignment2movie
        yes
        | ffmpeg
            -f image2
            -i snapshot%08d.png
            -pix_fmt yuv420p
            -vcodec libx264
            -profile:v baseline
            -level 30
            -maxrate 10000000
            -bufsize 10000000
            -b:v 1200k
            -f mp4
            movie.mp4)
