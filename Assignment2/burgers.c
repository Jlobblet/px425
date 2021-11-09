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

#ifdef __OPTIMIZE__
void allocate2d(double*** a, int Nx, int Ny) {
    double** b_loc;

    b_loc = (double**) calloc((Nx + 1) * Ny, sizeof(double*));
    if (b_loc == NULL) {
        fprintf(stderr, "malloc error in allocate2d\n");
        fflush(stderr);
        exit(1);
    }

    for (int iy = 0; iy < Nx; iy++) {

        b_loc[iy] = (double*) b_loc + (iy + 1) * Nx;
    }

    *a = b_loc;
}

void free2d(double*** a) {
    double** b_loc = *a;
    /* Release memory */
    free(b_loc);
    *a = NULL;
}
#else
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
#endif // __OPTIMIZE__

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
#ifdef NDEBUG
    int stepincr = 10;
#endif // NDEBUG
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
         * I interior */

        /* Corners */

        double u00 = u[0][0];
        double u10 = u[1][0];
        double uE0 = u[Nx - 1][0];
        double u01 = u[0][1];
        double u0E = u[0][Ny - 1];
        double ue0 = u[Nx - 2][0];
        double uE1 = u[Nx - 1][1];
        double uEE = u[Nx - 1][Ny - 1];
        double u1E = u[1][Ny - 1];
        double u0e = u[0][Ny - 2];

        /* Bottom left
         * ix = 0, iy = 0 */
        Lapl = (-2.0 * u00 + u10 + uE0) * inv_dx_sq /* left */
               + (-2.0 * u00 + u01 + u0E) * inv_dy_sq; /* bottom*/
        grad = (u10 - uE0) * half_inv_dx /* left */
               + (u01 - u0E) * half_inv_dy; /* bottom */
        u_new[0][0] = u00 - dt * u00 * grad + dt * nu * Lapl;

        /* Bottom right
         * ix = Nx - 1, iy = 0 */
        Lapl = (-2.0 * uE0 + u00 + ue0) * inv_dx_sq /* right */
               + (-2.0 * uE0 + uE1 + uEE) * inv_dy_sq;  /* bottom */
        grad = (u00 - ue0) * half_inv_dx /* right */
               + (uE1 - uEE) * half_inv_dy; /* bottom */
        u_new[Nx - 1][0] = uE0 - dt * uE0 * grad + dt * nu * Lapl;

        /* Top left
         * ix = 0, iy = Ny - 1 */
        Lapl = (-2.0 * u0E + u1E + uEE) * inv_dx_sq /* left */
               + (-2.0 * u0E + u00 + u0e) * inv_dy_sq; /* top */
        grad = (u1E - uEE) * half_inv_dx /* left */
               + (u00 - u0e) * half_inv_dx; /* top */
        u_new[0][Ny - 1] = u0E - dt * u0E * grad + dt * nu * Lapl;

        /* Top right
         * ix = Nx - 1, iy = Ny - 1 */
        Lapl = (-2.0 * uEE + u0E + u[Nx - 1 - 1][Ny - 1]) * inv_dx_sq /* right */
               + (-2.0 * uEE + uE0 + u[Nx - 1][Ny - 1 - 1]) * inv_dy_sq; /* top */
        grad = (u0E - u[Nx - 1 - 1][Ny - 1]) * half_inv_dx + /* right */
               +(uE0 - u[Nx - 1][Ny - 1 - 1]) * half_inv_dx; /* top */
        u_new[Nx - 1][Ny - 1] = uEE - dt * uEE * grad + dt * nu * Lapl;

        /* Edges */
        /* Left and right
         * ix = 0 and ix = Nx - 1 */
        for (iy = 1; iy < Ny - 1; iy++) {
            double u0y = u[0][iy];
            double u1y = u[1][iy];
            double uEy = u[Nx - 1][iy];
            double u0yp = u[0][iy + 1];
            double u0ym = u[0][iy - 1];
            double uey = u[Nx - 1 - 1][iy];
            double uEyp = u[Nx - 1][iy + 1];
            double uEym = u[Nx - 1][iy - 1];

            Lapl = (-2.0 * u0y + u1y + uEy) * inv_dx_sq  /* left */
                   + (-2.0 * u0y + u0yp + u0ym) * inv_dy_sq;
            grad = (u1y - uEy) * half_inv_dx
                   + (u0yp - u0ym) * half_inv_dy;
            u_new[0][iy] = u0y - dt * u0y * grad + dt * nu * Lapl;

            Lapl = (-2.0 * uEy + u0y + uey) * inv_dx_sq /* right */
                   + (-2.0 * uEy + uEyp + uEym) * inv_dy_sq;
            grad = (u0y - uey) * half_inv_dx
                   + (uEyp - uEym) * half_inv_dy;
            u_new[Nx - 1][iy] = uEy - dt * uEy * grad + dt * nu * Lapl;
        }
        /* Bottom and top
         * iy = 0 and iy = Ny - 1 */
        for (ix = 1; ix < Nx - 1; ix++) {
            double ux0 = u[ix][0];
            double uxp0 = u[ix + 1][0];
            double uxm0 = u[ix - 1][0];
            double ux1 = u[ix][0 + 1];
            double uxE = u[ix][Ny - 1];
            double uxpE = u[ix + 1][Ny - 1];
            double uxmE = u[ix - 1][Ny - 1];
            double uxe = u[ix][Ny - 1 - 1];

            Lapl = (-2.0 * ux0 + uxp0 + uxm0) * inv_dx_sq
                   + (-2.0 * ux0 + ux1 + uxE) * inv_dy_sq;  /* bottom */
            grad = (uxp0 - uxm0) * half_inv_dx
                   + (ux1 - uxE) * half_inv_dy;
            u_new[ix][0] = ux0 - dt * ux0 * grad + dt * nu * Lapl;

            Lapl = (-2.0 * uxE + uxpE + uxmE) * inv_dx_sq
                   + (-2.0 * uxE + ux0 + uxe) * inv_dy_sq;     /* top */
            grad = (uxpE - uxmE) * half_inv_dx
                   + (ux0 - uxe) * half_inv_dx;
            u_new[ix][Ny - 1] = uxE - dt * uxE * grad + dt * nu * Lapl;
        }

        /* Interior */
        for (ix = 1; ix < Nx - 1; ix++) {
            for (iy = 1; iy < Ny - 1; iy++) {
                double uxy = u[ix][iy];
                double uxpy = u[ix + 1][iy];
                double uxmy = u[ix - 1][iy];
                double uxyp = u[ix][iy + 1];
                double uxym = u[ix][iy - 1];

                Lapl = 0.0;
                grad = 0.0;

                /* Compute d2u/dx2 */
                Lapl += (-2.0 * uxy + uxpy + uxmy) * inv_dx_sq;
                grad += (uxpy - uxmy) * half_inv_dx;

                /* Compute d2u/dy2 */
                Lapl += (-2.0 * uxy + uxyp + uxym) * inv_dy_sq;
                grad += (uxyp - uxym) * half_inv_dy;

                /* Compute new value of u at this grid point */
                u_new[ix][iy] = uxy - dt * uxy * grad + dt * nu * Lapl;
            }
        }

        /* Shunt u_new into u */
        double** tmp = u_new;
        u_new = u;
        u = tmp;
#ifdef NDEBUG
        /* Snapshots of grid to file */
        if (istep == isnap) {
            sprintf(filename, "snapshot%08d.png", isnap);
            writePNG(filename, u, Nx, Ny);
            isnap *= stepincr;
        }
#endif // NDEBUG
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
#ifdef __OPTIMIZE__
    free2d(&u);
    free2d(&u_new);
#else
    free2d(&u, Nx);
    free2d(&u_new, Nx);
#endif // __OPTIMIZE__

    return EXIT_SUCCESS;
}
