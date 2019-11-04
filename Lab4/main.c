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
#include <gsl/gsl_linalg.h>
#include "dislin.h"

// Constant definitions
#define a -10.                        // lower bound
#define b  10.                        // upper bound
#define k  20                         // number of grid points
#define h (b-a)/(k-1)                 // distance between grid points
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

double solution(double *coeff, double x){
    assert(x >= a && x <= b);
    return 0;
}

// Main program
int main(int argc, char* argv[]) {
    printf("Program started with params"
           "\n[a,b] = [%f, %f]\n"
           "u(a) = %f\n"
           "u(b) = %f\n"
           "grid_point = %d\n", a, b, u_a, u_b, k);

    double test = c_3(1,1,1);
    printf("\n%f\n", test);
    return 0;
}
