/*
 * @author Dawid Majchrowski
 * @date   ?
 *
 * Galerkin method:
 * | u''(x) + p(x)*u'(x) + q(x)*u(x) = f(x)
 * | u(a) = u_a
 * | u(b) = u_b
*/

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include "dislin.h"

// Constant definitions
#define a -10.                        // lower bound
#define b  10.                        // upper bound
#define k  20                         // number of grid points - 2 Points: (0,1,...,k,k+1)
#define h (b-a)/(k+1)                 // distance between grid points
#define u_a  2.                       // u(a)
#define u_b  35.                      // u(b)
#define n 30                          // grid size for plotting

// Macros
#define MAX(A,B) (A) > (B) ? (A) : (B)
#define MIN(A,B) (A) > (B) ? (B) : (A)

// Method definitions
double p(double x);
double q(double x);
double f(double x);
double c_2 (int j, int i);
double cf_2(int j, int i);
double c_3 (int j, int m, int i);
double cf_3(int j, int m, int i);
double r(int j, int mode);
double *solve_linear_system(bool print_system, bool print_result);

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

double c_2(const int j, const int i){
    assert(abs(i-j) < 2);
    if(j==i)           return 2*h/3;
    return h/6;
}

double cf_2(const int j, const int i){
    assert(abs(i-j) < 2);
    if(j==i)           return -1/h;
    return  2/h;
}

double c_3(const int j, const int m, const int i){
    assert(abs(i-j) < 2 && abs(i-m) < 2 && abs(j-m) < 2);
    if(i==j && j==m) return h/2;
    return h/12;

}

double cf_3(const int j, const int m, const int i){
    assert(abs(i-j) < 2 && abs(i-m) < 2 && abs(j-m) < 2);
    if(j == m) {
        if(i==j)   return     0;
        if(m==i-1) return  1./3;
        if(m==i+1) return -1./3;
    }else if(j==i){
        if(i==m-1) return -1./6;
        if(i==m+1) return  1./6;
    }else {  // i==m
        if(i==j-1) return -1./6;
    }
    return  1./6;  // i==m  && i==j+1

}

double r(const int j, const int mode){  // compute value for sparse matrix coeff
    assert(j>=1 && j <= k);
    assert(mode >= -1 && mode <= 1); // r_j_j-1 or r_j_j or r_j_j+1
    assert(!(mode == -1 && j==1));
    assert(!(mode ==  1 && j==k));
    double x_j = a+j*h;
    double x_jm1 = x_j - h; // x_(j-1)
    double x_jp1 = x_j + h; // x_(j+1)

    if(mode == 1){       // r_j_j+1
            return -cf_2(j,j+1) + p(x_j)*cf_3(j,j,j+1) + q(x_j)*c_3(j,j,j+1) + p(x_jp1)*cf_3(j,j+1,j+1) + q(x_jp1)*c_3(j,j+1,j+1);
    }else if(mode == 0){ // r_j_j
            return -cf_2(j,j) + p(x_jm1)*cf_3(j,j-1,j) + q(x_jm1)*c_3(j,j-1,j) + p(x_j)*cf_3(j,j,j) + q(x_j)*c_3(j,j,j)
            + p(x_jp1)*cf_3(j,j+1,j) + q(x_jp1)*c_3(j,j+1,j);
        }else {          // r_j_j-1 && mode == -1
        return -cf_2(j, j - 1) + p(x_jm1) * cf_3(j, j - 1, j - 1) + q(x_jm1) * c_3(j, j - 1, j - 1) +
               p(x_j) * cf_3(j, j, j - 1) + q(x_j) * c_3(j, j, j - 1);
    }
}

double s(const int j){ // right side of equations
    assert(j >=1 && j<= k-1);
    double x_j = a+j*h;
    double x_jm1 = x_j - h; // x_(j-1)
    double x_jp1 = x_j + h; // x_(j+1)

    double result = f(x_jm1)*c_2(j,j-1) + f(x_j)*c_2(j,j) + f(x_jp1)*c_2(j,j+1);
    if(j == 1) result = result + u_a*(cf_2(1,0) - p(a)*cf_3(1,0,0) - p(a+h)*cf_3(1,1,0) - q(a)*c_3(1,0,0) - q(a+h)*c_3(1,1,0));
    if(j == k) result = result + u_b*(cf_2(k,k+1) - p(b-h)*cf_3(k,k,k+1) - p(b)*cf_3(k,k+1,k+1) - q(b-h)*c_3(k,k,k+1) - q(b)*c_3(k,k+1,k+1));
    return result;
}

//double *solve_linear_system(bool print_system, bool print_result){ // calculate alpha_i for i e {1,...,k}
//    double *coeffs = calloc(k, sizeof(double));
//    coeffs[0] = u_a;   // alpha(1)
//    coeffs[1] = u_b;   // alpha(2)
//
//    double *matrix = calloc((k-2)*(k-2), sizeof(double));  // (k-2)x(k-2) matrix of coefficients
//    double *values = calloc(k-2, sizeof(double));          //  k-2        array of function values
//
//    double x;
//    for(int i = 1; i < k-1; ++i){
//        x = a + i*h;
//        for(int deg = 2; deg < k; ++deg){
//            matrix[(i-1)*(k-2) + (deg-2)] = coeff(deg, x);
//        }
//        values[(i-1)] = f(x) - u_a * coeff(1, x) - u_b * coeff(2, x);
//        }
//
//    if(print_system){ // if flag is set, print the system of equations to be solved
//        printf("---------------- System of equations ----------------------------------\n");
//        for(int i = 0; i < k-2; ++i){
//            for(int j = 0; j < k-2; ++j){
//                printf("%10.6f ", matrix[i*(k-2)+j]);
//            }
//            printf("= %10.6f\n", values[i]);
//        }
//    }
//
//    gsl_matrix_view matrix_view = gsl_matrix_view_array (matrix, (k-2), (k-2));
//
//    gsl_vector_view vector_view = gsl_vector_view_array (values, k-2);
//
//    gsl_vector *result = gsl_vector_alloc (k-2);
//
//    int s;
//
//    gsl_permutation *p = gsl_permutation_alloc (k-2);
//
//    gsl_linalg_LU_decomp (&matrix_view.matrix, p, &s);
//
//    gsl_linalg_LU_solve (&matrix_view.matrix, p, &vector_view.vector, result);
//
//    for(int i = 0; i < k-2; ++i){
//        coeffs[i+2] = result->data[i];
//    }
//
//    if(print_result){ // if flag is set, print the solution
//        printf("---------------- Solution of given system of equations ----------------\n");
//        for(int i = 0; i < k; ++i){
//            printf("alpha_%d = %f\n", i+1, coeffs[i]);
//        }
//    }
//
//    gsl_permutation_free (p);
//    gsl_vector_free (result);
//    free(matrix);
//    free(values);
//    return coeffs;
//
//}

//void plot_result(const double *coeffs, const int N){
//    float *arguments = calloc(N, sizeof(float));
//    float *values    = calloc(N, sizeof(float));
//    float step = (b-a)/(N-1);
//    int ic =   intrgb (0.95,0.95,0.95);
//    double x_l = a - 5.*step;             // x_left
//    double x_r = b + 5.*step;             // x_right
//    double y_b = MIN(u_a, u_b) - 5.*step; // y_bottom
//    double y_t = MAX(u_a, u_b) + 5.*step; // y_top
//
//    for (int i = 0; i < N; ++i){
//        arguments[i] = a + i*step;
//        values[i] = final_u(arguments[i], coeffs);
//        y_b = MIN(y_b, values[i] - 5.*step); // update OY if needed
//        y_t = MAX(y_t, values[i] + 5.*step); // update OY if needed
//    }
//    metafl ("cons");
//    disini ();
//    pagera ();
//    complx ();
//    axspos (450, 1800);
//    axslen (2200, 1200);
//    name   ("OX", "x");
//    name   ("OY", "y");
//    labdig (1,  "xy");
//    ticks  (10, "x");
//    ticks  (2, "y");
//    titlin ("Collocation method solution", 1);
//    axsbgd (ic);
//    graf   (x_l, x_r, a, b/4., y_b, y_t, y_b, y_t/4.);
//    setrgb (0.7, 0.7, 0.7);
//    grid   (1, 1);
//    color  ("fore");
//    height (50);
//    title  ();
//    color  ("blue");
//    curve  (arguments, values, N);
//    disfin ();
//    free(arguments);
//    free(values);
//}


double solution(double *coeff, double x){
    assert(x >= a && x <= b);
    return 0;
}

 Main program
int main(int argc, char* argv[]) {
    printf("Program started with params"
           "\n[a,b] = [%f, %f]\n"
           "u(a) = %f\n"
           "u(b) = %f\n"
           "grid_point = %d\n", a, b, u_a, u_b, k);
    double *coeffs = solve_linear_system(true, true); // if true will print (system of equations/solution)
    plot_result(coeffs, n+1);
    return 0;
}
