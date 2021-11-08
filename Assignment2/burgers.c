/*  ANSI-C code (un optimised) for PX425 assignment 2 2021
 *  Evolves a function u in 2D via finite differences
 *  of the 2D Burgers Equation (simplified)
 *
 *  d u          du   du           d^2       d^2
 *  ----- = -u ( -- + -- ) + nu ( ---- u  +  ---- u )
 *  d t          dx   dy          dx^2       dy^2
 *
 *  Based originally on code created by D. Quigley
 *  Adapted by N. Hine */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "makePNG.h"
#include "mt19937ar.h"

/* Function prototypes for memory management routines */

void allocate2d(double*** a, int Nx, int Ny);

void free2d(double*** a, int Nx);

int main() {
    /* Function u on new and current grid */
    double** u_new, ** u;

    /* Approximate Laplacian */
    double Lapl, grad;

    /* Number of grid points */
    const int Nx = 256;
    const int Ny = 256;

    /* Loop counters */
    int ix, iy, istep;

    /* Filename to which the grid is drawn */
    int isnap = 0;
    char filename[25];

    /* Initial time */
    clock_t t1 = clock();

    /* Initialise random number generator */
    const unsigned long seed = 120549784972;
    init_genrand(seed);

    /* Set grid spacing */
    const double dx = 1.0;
    const double dy = 1.0;
    const double dx_sq = dx * dx;
    const double dy_sq = dy * dy;
    const double inv_dx = 1.0 / dx;
    const double inv_dy = 1.0 / dy;
    const double inv_dx_sq = 1.0 / dx_sq;
    const double inv_dy_sq = 1.0 / dy_sq;
    const double half_inv_dx = inv_dx * 0.5;
    const double half_inv_dy = inv_dy * 0.5;

    /* Set timestep and K */
    const double dt = 0.001;
    const double nu = 5.0;

    /* Number of steps to run */
    const int nstep = 10000;

    /* Allocate memory for a bunch of stuff */
    allocate2d(&u, Nx, Ny);
    allocate2d(&u_new, Nx, Ny);

    /* Initialise with random numbers */
    for (ix = 0; ix < Nx; ix++) {
        for (iy = 0; iy < Ny; iy++) {
            u[ix][iy] = 2.0 * genrand() - 1.0;
        }
    }

    /* Write an image of the initial grid */
    int stepincr = 10;
    sprintf(filename, "snapshot%08d.png", isnap);
    writePNG(filename, u, Nx, Ny);
    isnap++;

    /* Setup time */
    clock_t t2 = clock();
    printf("Setup time                    : %15.6f seconds\n", (double) (t2 - t1) / (double) CLOCKS_PER_SEC);
    fflush(stdout);
    t1 = t2;

    /* BEGIN SECTION TO BE OPTIMISED */

    /* Loop over the number of output timesteps */
    for (istep = 1; istep < nstep; istep++) {

        /* Quick diagram to label the upcoming sections:
         *
         * CEEEEC
         * EIIIIE
         * EIIIIE
         * EIIIIE
         * EIIIIE
         * CEEEEC
         *
         * C corner
         * E edge
         * I interior
         *
         * The compiler will take care of expressions like 0 + 1 in this bit. */

        /* Corners */

        /* Bottom left
         * ix = 0, iy = 0 */
        Lapl = (-2.0 * u[0][0] + u[0 + 1][0] + u[Nx - 1][0]) * inv_dx_sq /* left */
               + (-2.0 * u[0][0] + u[0][0 + 1] + u[0][Ny - 1]) * inv_dy_sq; /* bottom*/
        grad = (u[0 + 1][0] - u[Nx - 1][0]) * half_inv_dx /* left */
               + (u[0][0 + 1] - u[0][Ny - 1]) * half_inv_dy; /* bottom */
        u_new[0][0] = u[0][0] - dt * u[0][0] * grad + dt * nu * Lapl;

        /* Bottom right
         * ix = Nx - 1, iy = 0 */
        Lapl = (-2.0 * u[Nx - 1][0] + u[0][0] + u[Nx - 1 - 1][0]) * inv_dx_sq /* right */
               + (-2.0 * u[Nx - 1][0] + u[Nx - 1][0 + 1] + u[Nx - 1][Ny - 1]) * inv_dy_sq;  /* bottom */
        grad = (u[0][0] - u[Nx - 1 - 1][0]) * half_inv_dx /* right */
               + (u[Nx - 1][0 + 1] - u[Nx - 1][Ny - 1]) * half_inv_dy; /* bottom */
        u_new[Nx - 1][0] = u[Nx - 1][0] - dt * u[Nx - 1][0] * grad + dt * nu * Lapl;

        /* Top left
         * ix = 0, iy = Ny - 1 */
        Lapl = (-2.0 * u[0][Ny - 1] + u[0 + 1][Ny - 1] + u[Nx - 1][Ny - 1]) * inv_dx_sq /* left */
               + (-2.0 * u[0][Ny - 1] + u[0][0] + u[0][Ny - 1 - 1]) * inv_dy_sq; /* top */
        grad = (u[0 + 1][Ny - 1] - u[Nx - 1][Ny - 1]) * half_inv_dx /* left */
               + (u[0][0] - u[0][Ny - 1 - 1]) * half_inv_dx; /* top */
        u_new[0][Ny - 1] = u[0][Ny - 1] - dt * u[0][Ny - 1] * grad + dt * nu * Lapl;

        /* Top right
         * ix = Nx - 1, iy = Ny - 1 */
        Lapl = (-2.0 * u[Nx - 1][Ny - 1] + u[0][Ny - 1] + u[Nx - 1 - 1][Ny - 1]) * inv_dx_sq /* right */
               + (-2.0 * u[Nx - 1][Ny - 1] + u[Nx - 1][0] + u[Nx - 1][Ny - 1 - 1]) * inv_dy_sq; /* top */
        grad = (u[0][Ny - 1] - u[Nx - 1 - 1][Ny - 1]) * half_inv_dx + /* right */
               +(u[Nx - 1][0] - u[Nx - 1][Ny - 1 - 1]) * half_inv_dx; /* top */
        u_new[Nx - 1][Ny - 1] = u[Nx - 1][Ny - 1] - dt * u[Nx - 1][Ny - 1] * grad + dt * nu * Lapl;

        /* Edges */
        /* Left and right
         * ix = 0 and ix = Nx - 1 */
        for (iy = 1; iy < Ny - 1; iy++) {
            Lapl = (-2.0 * u[0][iy] + u[0 + 1][iy] + u[Nx - 1][iy]) * inv_dx_sq  /* left */
                   + (-2.0 * u[0][iy] + u[0][iy + 1] + u[0][iy - 1]) * inv_dy_sq;
            grad = (u[0 + 1][iy] - u[Nx - 1][iy]) * half_inv_dx
                   + (u[0][iy + 1] - u[0][iy - 1]) * half_inv_dy;
            u_new[0][iy] = u[0][iy] - dt * u[0][iy] * grad + dt * nu * Lapl;

            Lapl = (-2.0 * u[Nx - 1][iy] + u[0][iy] + u[Nx - 1 - 1][iy]) * inv_dx_sq /* right */
                   + (-2.0 * u[Nx - 1][iy] + u[Nx - 1][iy + 1] + u[Nx - 1][iy - 1]) * inv_dy_sq;
            grad = (u[0][iy] - u[Nx - 1 - 1][iy]) * half_inv_dx
                   + (u[Nx - 1][iy + 1] - u[Nx - 1][iy - 1]) * half_inv_dy;
            u_new[Nx - 1][iy] = u[Nx - 1][iy] - dt * u[Nx - 1][iy] * grad + dt * nu * Lapl;
        }
        /* Bottom and top
         * iy = 0 and iy = Ny - 1 */
        for (ix = 1; ix < Nx - 1; ix++) {
            Lapl = (-2.0 * u[ix][0] + u[ix + 1][0] + u[ix - 1][0]) * inv_dx_sq
                   + (-2.0 * u[ix][0] + u[ix][0 + 1] + u[ix][Ny - 1]) * inv_dy_sq;  /* bottom */
            grad = (u[ix + 1][0] - u[ix - 1][0]) * half_inv_dx
                   + (u[ix][0 + 1] - u[ix][Ny - 1]) * half_inv_dy;
            u_new[ix][0] = u[ix][0] - dt * u[ix][0] * grad + dt * nu * Lapl;

            Lapl = (-2.0 * u[ix][Ny - 1] + u[ix + 1][Ny - 1] + u[ix - 1][Ny - 1]) * inv_dx_sq
                   + (-2.0 * u[ix][Ny - 1] + u[ix][0] + u[ix][Ny - 1 - 1]) * inv_dy_sq;     /* top */
            grad = (u[ix + 1][Ny - 1] - u[ix - 1][Ny - 1]) * half_inv_dx
                   + (u[ix][0] - u[ix][Ny - 1 - 1]) * half_inv_dx;
            u_new[ix][Ny - 1] = u[ix][Ny - 1] - dt * u[ix][Ny - 1] * grad + dt * nu * Lapl;
        }

        /* Interior */
        for (ix = 1; ix < Nx - 1; ix++) {
            for (iy = 1; iy < Ny - 1; iy++) {
                Lapl = 0.0;
                grad = 0.0;

                /* Compute d2u/dx2 */
                Lapl += (-2.0 * u[ix][iy] + u[ix + 1][iy] + u[ix - 1][iy]) * inv_dx_sq;
                grad += (u[ix + 1][iy] - u[ix - 1][iy]) * half_inv_dx;

                /* Compute d2u/dy2 */
                Lapl += (-2.0 * u[ix][iy] + u[ix][iy + 1] + u[ix][iy - 1]) * inv_dy_sq;
                grad += (u[ix][iy + 1] - u[ix][iy - 1]) * half_inv_dy;

                /* Compute new value of u at this grid point */
                u_new[ix][iy] = u[ix][iy] - dt * u[ix][iy] * grad + dt * nu * Lapl;
            }
        }

        /* Shunt u_new into u */
        double** tmp = u_new;
        u_new = u;
        u = tmp;

        /* Snapshots of grid to file */
        if (istep == isnap) {
            sprintf(filename, "snapshot%08d.png", isnap);
            writePNG(filename, u, Nx, Ny);
            isnap *= stepincr;
        }

    }

    /* END SECTION TO BE OPTIMISED */

    /* Calculation time */
    t2 = clock();
    printf("Time taken for %8d steps : %15.6f seconds\n", nstep, (double) (t2 - t1) / (double) CLOCKS_PER_SEC);
    fflush(stdout);

    /* Write an image of the final grid */
    sprintf(filename, "snapshot%08d.png", istep);
    writePNG(filename, u, Nx, Ny);

    /* Write final time-evolved solution to file */
    FILE* fp = fopen("final_grid.dat", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening final_grid.dat for output\n");
        fflush(stderr);
    }

    for (ix = 0; ix < Nx - 1; ix++) {
        for (iy = 0; iy < Ny - 1; iy++) {
            /* x and y at the current grid points */
            double x = dx * (double) ix;
            double y = dy * (double) iy;
            fprintf(fp, "%8.4f %8.4f %8.4e\n", x, y, u[ix][iy]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    /* Release memory */
    free2d(&u, Nx);
    free2d(&u_new, Nx);

    return EXIT_SUCCESS;
}


/* Auxiliary routines for memory management */

void allocate2d(double*** a, int Nx, int Ny) {

    double** b_loc;

    b_loc = (double**) calloc(Nx, sizeof(double*));
    if (b_loc == NULL) {
        fprintf(stderr, "malloc error in allocate2d\n");
        fflush(stderr);
    }

    int iy;
    for (iy = 0; iy < Nx; iy++) {

        b_loc[iy] = (double*) calloc(Ny, sizeof(double));
        if (b_loc[iy] == NULL) {
            fprintf(stderr, "malloc error for row %d of %d in allocate2d\n", iy, Nx);
            fflush(stderr);
        }

    }

    *a = b_loc;
}

void free2d(double*** a, int Nx) {
    int iy;
    double** b_loc = *a;
    /* Release memory */
    for (iy = 0; iy < Nx; iy++) {
        free(b_loc[iy]);
    }
    free(b_loc);
    *a = NULL;
}
