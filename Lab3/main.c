/*
 * @author Dawid Majchrowski
 * @date   22.10.2019
 *
 * Collocation method for given system of equations:
 * | u''(x) + p(x)*u'(x) + q(x)*u(x) = f(x)
 * | u(a) = u_a
 * | u(b) = u_b
 *
 * For given example
 * [a,b] = [-1, 2]
 * u(x) = x^5 + 3, p(x) = x^2 * cos(x), q(x) = x * sin(x)
 * Other values are calculated analytically
*/

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_linalg.h>
#include "dislin.h"

// Constant definitions
#define a -1.                         // lower bound
#define b  2.                         // upper bound
#define k  4                          // number of grid points
#define h (b-a)/(k-1)                 // distance between grid points
#define u_a  2.                       // u(a)
#define u_b  35.                      // u(b)
#define derivative "analytical"       // "empirical" or "analytical" calculation of Peano derivative
#define n 30                          // grid size for plotting

// Macros
#define MAX(A,B) (A) > (B) ? (A) : (B)
#define MIN(A,B) (A) > (B) ? (B) : (A)

// Method definitions
double p(double x);
double q(double x);
double f(double x);
double peano_polynomial(int deg, double x);
double analytical_first_derivative(int deg, double x);
double analytical_second_derivative(int deg, double x);
double empirical_first_derivative(int deg, double x);
double empirical_second_derivative(int deg, double x);
double peano_first_derivative(int deg, double x, const char* mode);
double peano_second_derivative(int deg, double x, const char* mode);
double coeff(int deg, double x);
double final_u(double x, const double *coeffs);
void plot_result(const double *coeffs, int N);

// Method implementations
double p(const double x){
    return x*x*cos(x);
}

double q(const double x){
    return x*sin(x);
}

double f(const double x){
    return 20*pow(x, 3) + 5*p(x)*pow(x, 4) + q(x)*(pow(x, 5) + 3);
}

double peano_polynomial(const int deg, const double x){
    assert(x > a - 10e-4  && x < b + 10e-4);   // t_k e [a, ... ,b]
    assert(deg >= 1 && deg <=k);               // deg e {1, 2, .. k}

    switch (deg) {
        case 1:  return (b - x)/(b - a);
        case 2:  return (x - a)/(b - a);
        case 3:  return ((b - x)*(x - a))/((b - a)*(b - a));
        default: return ((b - x)*(x - a)*pow((2*x - a - b), deg - 3))/pow((b - a), deg - 1);
    }
}

double analytical_first_derivative(const int deg, const double x){
    switch (deg) {
        case 1:  return (-1.)/(b - a);
        case 2:  return (1.) /(b - a);
        case 3:  return (b + a - 2*x)/((b - a)*(b - a));
        default: return ((pow((2*x - a - b), deg - 2))*((2*deg - 6)*(b - x)*(x - a) - (2*x - a - b)*(2*x - a - b)))/pow((b - a), deg - 1);
    }
}

double analytical_second_derivative(const int deg, const double x){
    switch (deg) {
        case 1:  return 0;
        case 2:  return 0;
        case 3:  return (-2.)/((b - a)*(b - a));
        case 4:  return (6.*(a + b - 2*x))/((b - a)*(b - a)*(b - a));
        default: return ((10 - 4*deg)*(pow((2*x - a - b), deg - 3)) + 4*(deg - 3)*(deg - 4)*(b - x)*(x - a)*pow((2*x - a - b), deg - 5))/pow((b - a), deg - 1);
    }
}

double empirical_first_derivative(const int deg, const double x){
    double eps = MIN(h, 10e-4);
    double next = x + eps;
    double prev = x - eps;
    return (peano_polynomial(deg, next) - peano_polynomial(deg, prev))/(2*eps);
}

double empirical_second_derivative(const int deg, const double x){
    double eps = MIN(h, 10e-4);
    double next = x + eps;
    double prev = x - eps;
    return (peano_polynomial(deg, next) + peano_polynomial(deg, prev) - 2*peano_polynomial(deg, x))/(eps*eps);
}

double peano_first_derivative(const int deg, const double x, const char* mode){
    assert(!strcmp(mode, "analytical") || ! strcmp(mode, "empirical")); // mode e {"analytical", "empirical"}
    assert(x > a && x < b);                                             // t_k  e [a + h, ... , b - h]
    assert(deg >= 1 && deg <=k);                                        // deg  e {1, 2, .. k}

    double (*first_derivative)(const int, const double) = !strcmp(mode, "analytical") ? &analytical_first_derivative : &empirical_first_derivative;
    return first_derivative(deg, x);
}

double peano_second_derivative(const int deg, const double x, const char* mode){
    assert(!strcmp(mode, "analytical") || !strcmp(mode, "empirical"));  // mode e {"analytical", "empirical"}
    assert(x > a && x < b);                                             // t_k  e [a + h, ... , b - h]
    assert(deg >= 1 && deg <=k);                                        // deg  e {1, 2, .. k}

    double (*second_derivative)(const int, const double) = !strcmp(mode, "analytical") ? &analytical_second_derivative : &empirical_second_derivative;
    return second_derivative(deg, x);
}

double coeff(const int deg, const double x){ // calculate one coeff in system of equations for given grid u
    return        peano_second_derivative(deg, x, derivative) +
           p(x) * peano_first_derivative (deg, x, derivative) +
           q(x) * peano_polynomial(deg, x);
}

double *solve_linear_system(bool print_system, bool print_result){ // calculate alpha_i for i e {1,...,k}
    double *coeffs = calloc(k, sizeof(double));
    coeffs[0] = u_a;   // alpha(1)
    coeffs[1] = u_b;   // alpha(2)

    double *matrix = calloc((k-2)*(k-2), sizeof(double));  // (k-2)x(k-2) matrix of coefficients
    double *values = calloc(k-2, sizeof(double));          //  k-2        array of function values

    double x;
    for(int i = 1; i < k-1; ++i){
        x = a + i*h;
        for(int deg = 2; deg < k; ++deg){
            matrix[(i-1)*(k-2) + (deg-2)] = coeff(deg, x);
        }
        values[(i-1)] = f(x) - u_a * coeff(1, x) - u_b * coeff(2, x);
        }

    if(print_system){ // if flag is set, print the system of equations to be solved
        printf("---------------- System of equations ----------------------------------\n");
        for(int i = 0; i < k-2; ++i){
            for(int j = 0; j < k-2; ++j){
                printf("%10.6f ", matrix[i*(k-2)+j]);
            }
            printf("= %10.6f\n", values[i]);
        }
    }

    gsl_matrix_view matrix_view = gsl_matrix_view_array (matrix, (k-2), (k-2));

    gsl_vector_view vector_view = gsl_vector_view_array (values, k-2);

    gsl_vector *result = gsl_vector_alloc (k-2);

    int s;

    gsl_permutation *p = gsl_permutation_alloc (k-2);

    gsl_linalg_LU_decomp (&matrix_view.matrix, p, &s);

    gsl_linalg_LU_solve (&matrix_view.matrix, p, &vector_view.vector, result);

    for(int i = 0; i < k-2; ++i){
        coeffs[i+2] = result->data[i];
    }

    if(print_result){ // if flag is set, print the solution
        printf("---------------- Solution of given system of equations ----------------\n");
        for(int i = 0; i < k; ++i){
            printf("alpha_%d = %f\n", i+1, coeffs[i]);
        }
    }

    gsl_permutation_free (p);
    gsl_vector_free (result);
    free(matrix);
    free(values);
    return coeffs;

}

double final_u(const double x, const double *coeffs){
    double sum = 0;
    for(int i = 0; i < k; ++i){
        sum = sum + coeffs[i]*peano_polynomial(i+1, x); // coeffs[i] stays for alpha(i+1)
    }
    return sum;
}

void plot_result(const double *coeffs, const int N){
    float *arguments = calloc(N, sizeof(float));
    float *values    = calloc(N, sizeof(float));
    float step = (b-a)/(N-1);
    int ic =   intrgb (0.95,0.95,0.95);
    double x_l = a - 5.*step;             // x_left
    double x_r = b + 5.*step;             // x_right
    double y_b = MIN(u_a, u_b) - 5.*step; // y_bottom
    double y_t = MAX(u_a, u_b) + 5.*step; // y_top

    for (int i = 0; i < N; ++i){
        arguments[i] = a + i*step;
        values[i] = final_u(arguments[i], coeffs);
        y_b = MIN(y_b, values[i] - 5.*step); // update OY if needed
        y_t = MAX(y_t, values[i] + 5.*step); // update OY if needed
    }
    metafl ("cons");
    disini ();
    pagera ();
    complx ();
    axspos (450, 1800);
    axslen (2200, 1200);
    name   ("OX", "x");
    name   ("OY", "y");
    labdig (1,  "xy");
    ticks  (10, "x");
    ticks  (2, "y");
    titlin ("Collocation method solution", 1);
    axsbgd (ic);
    graf   (x_l, x_r, a, b/4., y_b, y_t, y_b, y_t/4.);
    setrgb (0.7, 0.7, 0.7);
    grid   (1, 1);
    color  ("fore");
    height (50);
    title  ();
    color  ("blue");
    curve  (arguments, values, N);
    disfin ();
    free(arguments);
    free(values);
}

// Main program
int main(int argc, char* argv[]) {
    printf("Program started with params"
           "\n[a,b] = [%f, %f]\n"
           "u(a) = %f\n"
           "u(b) = %f\n"
           "grid_point = %d\n", a, b, u_a, u_b, k);
    double *coeffs = solve_linear_system(true, true); // if true will print (system of equations/solution)
    plot_result(coeffs, n+1);
    free(coeffs);
    return 0;
}
