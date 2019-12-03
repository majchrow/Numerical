/*
 * @author Dawid Majchrowski
 * @date   3.12.2019
 *
 * degenerate kernel method for Fredholm integral equations
 * | L*u(x) - inf_{a}^{b}K(x,s)u(s)ds = g(x)
 * | Expected solution for given problem
 * | a,b, L, N, kij - given arbitrary
 * | g(x) calculated for given u(x) = x^2
 * | u(x) = x^2 (Expected result)
*/

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include "dislin.h"

// Constant definitions
#define a   -1.                       // lower bound
#define b    1.                       // upper bound
#define L    1.                       // constant lambda in given equation
#define N    5                        // number of grid points - 2 Points: (0,1,...,k,k+1)
#define n    1000                     // grid size for plotting
#define kij  1.
#define mode "approximation"          // approximation or analytic

// Macros
#define MAX(A, B) (A) > (B) ? (A) : (B)
#define MIN(A, B) (A) > (B) ? (B) : (A)

typedef struct {
    int i;
    int k;
} points;

// Method definitions
double alpha(int i, double x);
double beta(int i, double x);
double int_g(double x);
double g(double x);
double K(double x, double s);
double m(int k, int j);
double p(int k);
double a_mul_b(double x, void *params);
double g_mul_b(double x, void *params);
double int_ab(int i, int k);
double int_gb(int k);
double final_u(double x, const double *coeffs);
double *solve_linear_system(bool print_system, bool print_result);
void plot_result(const double *coeffs, int k);

// Method implementations
double alpha(const int i, const double x) {
    return pow(x, i);
}

double beta(const int i, const double x) {
    return pow(x, i);
}

double g(const double x) {
    if (!strcmp(mode, "approximation")) {
        return x * x - int_g(x);
    } else { // analytic calculation of g (using horner, coefs calculated analytically)
        double sum = 142. / 105;
        for (int i = 1; i < N; ++i) {
            sum = sum * x + (142. / 105);
        }
        return x * x - sum * x;
    }
}

double int_ku(double s, void *params) {
    double x = *(double *) params;
    return K(x, s) * s * s;
}

double int_g(const double x) {
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
    double result, error;
    double epsabs = 1e-7;
    double epsrel = 1e-7;
    size_t limit = 1000;
    gsl_function F;
    F.function = &int_ku;
    double params = x;
    F.params = &params;
    gsl_integration_qags(&F, a, b, epsabs, epsrel, limit, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);
    return result;
}

double K(const double x, const double s) {
    double sum = 0.;
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            sum = sum + kij * alpha(i, x) * beta(j, s);
        }
    }
    return sum;
}

double a_mul_b(double x, void *params) {
    points point = *(points *) params;
    return alpha(point.i, x) * beta(point.k, x);
}

double g_mul_b(double x, void *params) {
    int k = *(int *) params;
    return g(x) * beta(k, x);
}

double int_ab(int i, int k) {
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
    double result, error;
    double epsabs = 1e-7;
    double epsrel = 1e-7;
    size_t limit = 1000;
    gsl_function F;
    F.function = &a_mul_b;
    points params = {.i = i, .k = k};
    F.params = &params;
    gsl_integration_qags(&F, a, b, epsabs, epsrel, limit, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);
    return result;
}

double int_gb(int k) {
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
    double result, error;
    double epsabs = 1e-7;
    double epsrel = 1e-7;
    size_t limit = 1000;
    gsl_function F;
    F.function = &g_mul_b;
    int params = k;
    F.params = &params;
    gsl_integration_qags(&F, a, b, epsabs, epsrel, limit, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);
    return result;
}

double m(const int k, const int j) {
    double result = 0.;
    for (int i = 1; i <= N; ++i) {
        result = result + kij * int_ab(i, k);
    }
    return result;
}

double p(const int k) {
    return int_gb(k);
}

double *solve_linear_system(bool print_system, bool print_result) { // calculate c_i for i e {1,...,k}
    double *coeffs = calloc(N, sizeof(double));

    double *matrix = calloc(N * N, sizeof(double *));    // (N)x(N) matrix of coefficients
    double *values = calloc(N, sizeof(double));          //  N        array of function values

    for (int k = 1; k <= N; ++k) {
        for (int j = 1; j <= N; ++j) {
            matrix[(k - 1) * N + j - 1] = k == j ? 1 - m(k, j) : -m(k, j);
        }
        values[k - 1] = p(k);
    }

    if (print_system) { // if flag is set, print the system of equations to be solved
        printf("---------------- System of equations ----------------------------------\n");
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                printf("%10.6f ", matrix[i * N + j]);
            }
            printf("= %10.6f\n", values[i]);
        }
    }

    gsl_matrix_view matrix_view = gsl_matrix_view_array(matrix, N, N);

    gsl_vector_view vector_view = gsl_vector_view_array(values, N);

    gsl_vector *result = gsl_vector_alloc(N);

    int s;

    gsl_permutation *p = gsl_permutation_alloc(N);

    gsl_linalg_LU_decomp(&matrix_view.matrix, p, &s);

    gsl_linalg_LU_solve(&matrix_view.matrix, p, &vector_view.vector, result);

    for (int i = 0; i < N; ++i) {
        coeffs[i] = result->data[i];
    }

    if (print_result) { // if flag is set, print the solution
        printf("---------------- Solution of given system of equations ----------------\n");
        for (int i = 0; i < N; ++i) {
            printf("c_%d = %f\n", i + 1, coeffs[i]);
        }
    }

    gsl_permutation_free(p);
    gsl_vector_free(result);
    free(matrix);
    free(values);
    return coeffs;
}

double final_u(const double x, const double *coeffs) {
    double u = 0.;
    for (int i = 1; i <= N; ++i) {
        double inner = 0.;
        for (int j = 1; j <= N; ++j) {
            inner = inner + kij * coeffs[j - 1];
        }
        u = u + alpha(i, x) * inner;
    }
    return (1 / L) * (u + g(x));
}

void plot_result(const double *coeffs, const int k) {
    float *arguments = calloc(k, sizeof(double));
    float *values = calloc(k, sizeof(double));
    double step = (b - a) / (k - 1);
    int ic = intrgb(0.95, 0.95, 0.95);
    double x_l = a - 5. * step;             // x_left
    double x_r = b + 5. * step;             // x_right
    double y_b = -5. * step; // y_bottom
    double y_t = 5. * step; // y_top

    for (int i = 0; i < k; ++i) {
        arguments[i] = a + i * step;
        values[i] = final_u(arguments[i], coeffs);
        y_b = MIN(y_b, values[i] - 5. * step); // update OY if needed
        y_t = MAX(y_t, values[i] + 5. * step); // update OY if needed
    }
    metafl("cons");
    disini();
    pagera();
    complx();
    axspos(450, 1800);
    axslen(2200, 1200);
    name("OX", "x");
    name("OY", "y");
    labdig(1, "xy");
    ticks(10, "x");
    ticks(2, "y");
    titlin("?? method solution", 1);
    axsbgd(ic);
    graf(x_l, x_r, a, b / 4., y_b, y_t, y_b, y_t / 4.);
    setrgb(0.7, 0.7, 0.7);
    grid(1, 1);
    color("fore");
    height(50);
    title();
    color("blue");
    curve(arguments, values, k);
    disfin();
    free(arguments);
    free(values);
}

// Main program
int main(int argc, char *argv[]) {
    printf("Program started with params"
           "\n[a,b] = [%f, %f]\n"
           "N = %d\n", a, b, N);
    double *coeffs = solve_linear_system(true, true); // if true will print (system of equations/solution)
    plot_result(coeffs, n + 1);
    return 0;
}
