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
    int Nx = 256;
    int Ny = 256;

    /* Loop counters */
    int ix, iy, istep;

    /* Filename to which the grid is drawn */
    int isnap = 0;
    char filename[25];

    /* Initial time */
    clock_t t1 = clock();

    /* Initialise random number generator */
    unsigned long seed = 120549784972;
    init_genrand(seed);

    /* Set grid spacing */
    double dx = 1.0;
    double dy = 1.0;
    double dxsq = dx * dx;
    double dysq = dy * dy;
    double invdx = 1.0 / dx;
    double invdy = 1.0 / dy;
    double invdxsq = 1.0 / dxsq;
    double invdysq = 1.0 / dysq;

    /* Set timestep and K */
    double dt = 0.001;
    double nu = 5.0;

    /* Number of steps to run */
    int nstep = 10000;

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

        /* Loop over grid points */
        for (ix = 0; ix < Nx; ix++) {
            for (iy = 0; iy < Ny; iy++) {

                /* Initialise finite-difference Laplacian */
                Lapl = 0.0;
                /* Initialise gradient */
                grad = 0.0;

                /* Compute d2u/dx2, accounting for special cases near boundaries */
                if (ix == 0) {
                    Lapl += (-2.0 * u[ix][iy] + u[ix + 1][iy] + u[Nx - 1][iy]) * invdxsq;  /* left */
                    grad += (u[ix + 1][iy] - u[Nx - 1][iy]) * invdx * 0.5;
                } else if (ix == Nx - 1) {
                    Lapl += (-2.0 * u[ix][iy] + u[0][iy] + u[ix - 1][iy]) / dxsq;     /* right */
                    grad += (u[0][iy] - u[ix - 1][iy]) * invdx * 0.5;
                } else {
                    Lapl += (-2.0 * u[ix][iy] + u[ix + 1][iy] + u[ix - 1][iy]) / dxsq;
                    grad += (u[ix + 1][iy] - u[ix - 1][iy]) * invdx * 0.5;
                }

                /* Compute d2u/dy2, accounting for special cases near boundaries */
                if (iy == 0) {
                    Lapl += (-2.0 * u[ix][iy] + u[ix][iy + 1] + u[ix][Ny - 1]) * invdysq;  /* bottom */
                    grad += (u[ix][iy + 1] - u[ix][Ny - 1]) * invdy * 0.5;
                } else if (iy == Ny - 1) {
                    Lapl += (-2.0 * u[ix][iy] + u[ix][0] + u[ix][iy - 1]) * invdysq;     /* top */
                    grad += (u[ix][0] - u[ix][iy - 1]) * invdy * 0.5;
                } else {
                    Lapl += (-2.0 * u[ix][iy] + u[ix][iy + 1] + u[ix][iy - 1]) * invdysq;
                    grad += (u[ix][iy + 1] - u[ix][iy - 1]) * invdy * 0.5;
                }

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
