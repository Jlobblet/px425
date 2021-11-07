// Compiled and run on godzilla successfully.
// Compilation instructions:
// gcc -Wall -Wpedantic -Wextra -Werror -std=c99 -lm -o qho qho.c
// Example usage:
// ./qho <<< "0.0001 1000"
// Example output:
// Numerical integral:	0.000138155332545
// Analytic integral:	0.000138155332545
// Difference (error):	0.000000000000000

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>

// Pi is only defined in GNU extensions of math.h, so it is defined here
const double PI = 3.14159265358979323846;

typedef struct {
    double a;
    int i_max;
} input_params_t;

double analytic_result();

double calc_result(input_params_t* params);

double simpson(input_params_t* params);

double calc_f(double x);

double calc_x(int i, double a);

double calc_h(int i, double a);

bool read_input(input_params_t* params);

int main() {
    input_params_t params;
    if (!read_input(&params)) {
        return EXIT_FAILURE;
    }

    double r1 = calc_result(&params);
    double r2 = analytic_result();
    double difference = fabs(r1 - r2);
    printf("Numerical integral:\t%.15f\nAnalytic integral:\t%.15f\nDifference (error):\t%.15f\n", r1, r2, difference);

    return EXIT_SUCCESS;
}

/// The analytical solution to the integral
double analytic_result() {
    double mod_A_3 = 48 * sqrt(PI);
    return 1 / (mod_A_3 * mod_A_3);
}

/// Calculate 1/|A_3|^2
double calc_result(input_params_t* params) {
    double mod_A_3 = fabs(simpson(params));
    return 1 / (mod_A_3 * mod_A_3);
}

/// Implementation of Composite Simpson's Rule for irregularly spaced data
double simpson(input_params_t* params) {
    // Using Kahan summation here, aka compensated summation.
    // This is an algorithm for reducing error when summing finite-precision floating point numbers.
    // Since Kahan summation is designed to reduce error in the least significant bits, it may be overkill here.
    double sum = 0, c = 0, y;
    // t and z are marked as volatile to prevent compiler optimisations from destroying the algorithm.
    volatile double t, z;
    for (int i = 0; i < params->i_max / 2 - 1; i++) {
        // Building blocks for the summation term

        /// 2i
        int two_i = i * 2;
        /// x_{2i}
        double x_2i = calc_x(two_i, params->a);
        /// x_{2i+1}
        double x_2ip = calc_x(two_i + 1, params->a);
        /// x_{2i+2}
        double x_2ipp = calc_x(two_i + 2, params->a);
        /// h_{2i}
        double h_2_i = calc_h(two_i, params->a);
        /// h_{2i+1}
        double h_2_ip = calc_h(two_i + 1, params->a);
        // Summation terms
        /// Outer term
        /// (h_{2i} + h_{2i+1}) / 6
        double outer = (h_2_i + h_2_ip) / 6;
        /// Term 1
        /// (2 - h_{2i+1}/h_{2i}) f(x_{2i})
        double t1 = (2 - h_2_ip / h_2_i) * calc_f(x_2i);

        // Term 2

        /// (h_{2i} + h{2i+1})
        double hsum = h_2_i + h_2_ip;
        /// (h_{2i} + h{2i+1})^2
        double hsum_sq = hsum * hsum;
        /// (h_{2i} _ h_{2i+1})^2/h_{2i}h_{2i+1} f(x_{2i+1})
        double t2 = hsum_sq / (h_2_i * h_2_ip) * calc_f(x_2ip);

        // Term 3
        /// (2 - h_{2i} / h_{2i+1}) f(x_{2i+2})
        double t3 = (2 - h_2_i / h_2_ip) * calc_f(x_2ipp);

        // This is the overall contribution of this term of summation.
        double contribution = outer * (t1 + t2 + t3);

        // If the contribution is sufficiently small to not affect the output
        // (15 decimal places) then stop looping
        if (contribution < 1e-16) {
            break;
        }

        // Now perform Kahan summation with this and the running sum so far.
        y = contribution - c;
        t = sum + y;
        z = t - sum;
        c = z - y;
        sum = t;
    }
    return sum;
}

/// H_3(calc_x) is the Hermite polynomial of order 3, 8x^3 - 12x
double H_3(double x) {
    return x * (8 * x * x - 12);
}

/// calc_f(x) is the integrand of the expression, |H_3(x)|^2 exp(-calc_x^2)
double calc_f(double x) {
    /// H_3(x)
    double H_3_x = H_3(x);
    /// x^2
    double xsq = x * x;
    /// |H_3(x)|
    double mod_H_3_x = fabs(H_3_x);
    // Move the 2 inside the integral now, so I don't forget it later
    // Negligible performance impact anyway.
    /// 2 |H_3(x)|^2 exp(-x^2)
    return 2 * mod_H_3_x * mod_H_3_x * exp(-xsq);
}

/// x_i = ai^2;
double calc_x(int i, double a) {
    return a * i * i;
}

/// h_i = x_{i+1} - x_i = a(2i + 1)
double calc_h(int i, double a) {
    return a * (2 * i + 1);
}

/// Read input data consisting of a double-precision float and a long
bool read_input(input_params_t* params) {
    bool cont = true;
    char* start_ptr, * end_ptr;
    // 256 bytes should be enough to read in both numbers on one line
    char buf[256];
    // Using getline here for greater control over the input stream.
    if (fgets(buf, 256, stdin) == NULL) {
        fprintf(stderr, "Failed to read input.\n");
        cont = false;
    }
    // Using strtod and strtol so that if they fail it can be detected immediately.
    if (cont) {
        start_ptr = buf;
        params->a = strtod(start_ptr, &end_ptr);
        if (start_ptr == end_ptr) {
            // start = end means that the function failed to advance the stream at all, so nothing was parsed.
            fprintf(stderr, "Failed to parse a.\n");
            cont = false;
        }
    }
    // If we have reached the end of the line, read a new one in
    if (cont) {
        // Skip trailing whitespace
        while (isspace(*end_ptr)) {
            end_ptr++;
        }
        // End of the line
        if (*end_ptr == 0) {
            // Overwrite the existing buffer
            if (fgets(buf, 256, stdin) == NULL) {
                fprintf(stderr, "Failed to read input.\n");
                cont = false;
            } else {
                // Reset the pointer to the beginning
                end_ptr = buf;
            }
        }
    }
    // Now read in i_max
    if (cont) {
        // The start pointer here is where the previous parse attempt stopped.
        start_ptr = end_ptr;
        // strtol makes you specify the base (10).
        params->i_max = (int) strtol(start_ptr, &end_ptr, 10);
        if (start_ptr == end_ptr) {
            // Same as above.
            fprintf(stderr, "Failed to parse i_max.\n");
            cont = false;
        }
    }
    // If cont = false above, the code will short-circuit to here and return false.
    // If everything succeeds, the code will return true.
    // This style of programming allows for minimal handling around earlier return statements, e.g.
    // line only needs to be freed here and not in several places.
    return cont;
}
