/* Metoda kolokacji dla układu
 * | u''(x) + p(x)*u'(x) + q(x)*u = f(final_u)
 * | u(a) = u_a
 * | u(b) = u_b
 * Dawid Majchrowski 22.10.2019
*/

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

// Constant definitions
#define a 0                    // lower bound
#define b 19                   // upper bound
#define k 19                   // number of grid points
#define h 1                    // (b-a)/k
#define u_a 0                  // u(a)
#define u_b 0                  // u(b)
#define derivative "empirical" // empirical or analytical calculation of peano derivative

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
double coeff(double x, int deg);
double final_u(double x, const double *coeffs);

// Method implementations
double p(const double x){
    return x*x + 1;
}

double q(const double x){
    return x*x + cos(x);
}

double f(const double x){
    return x*x*sin(x);
}

double peano_polynomial(const int deg, const double x){
    assert(x >= a + h && x <= b - h);  // t_k    e [a + h, ... , b - h]
    assert(deg >= 1 && deg <=k);       // deg e {1, 2, .. k}
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
    double next = x + h;
    double prev = x - h;
    return (peano_polynomial(deg, next) - peano_polynomial(deg, prev))/(2*h);
}

double empirical_second_derivative(const int deg, const double x){
    double next = x + h;
    double prev = x - h;
    return (peano_polynomial(deg, next) + peano_polynomial(deg, prev) - 2*peano_polynomial(deg, x))/(h*h);
}

double peano_first_derivative(const int deg, const double x, const char* mode){
    assert(!strcmp(mode, "analytical") || ! strcmp(mode, "empirical")); // mode e {"analytical", "empirical"}
    assert(x >= a + h && x <= b - h);                                   // t_k    e [a + h, ... , b - h]
    assert(deg >= 1 && deg <=k);                                        // deg e {1, 2, .. k}

    double (*first_derivative)(const int, const double) = !strcmp(mode, "analytical") ? &analytical_first_derivative : &empirical_first_derivative;
    return first_derivative(deg, x);
}

double peano_second_derivative(const int deg, const double x, const char* mode){
    assert(!strcmp(mode, "analytical") || ! strcmp(mode, "empirical")); // mode e {"analytical", "empirical"}
    assert(x >= a + h && x <= b - h);                                   // t_k    e [a + h, ... , b - h]
    assert(deg >= 1 && deg <=k);                                        // deg e {1, 2, .. k}

    double (*second_derivative)(const int, const double) = !strcmp(mode, "analytical") ? &analytical_second_derivative : &empirical_second_derivative;
    return second_derivative(deg, x);
}

double coeff(const double x, const int deg){ // calculate one coeff in system of equations for given grid final_u
    return        peano_second_derivative(deg, x, derivative) +
           p(x) * peano_first_derivative (deg, x, derivative) +
           q(x) * peano_polynomial(deg, x);
}

double *solve_linear_system(){ // calculate alpha_i for i e {1,...,k}
    double *coeffs = calloc(k, sizeof(double));
    coeffs[0] = u_a; // alpha(1)
    coeffs[k-1] = u_b; // alpha(k)

    double a_data[(k-2)*(k-2)];
    double b_data[(k-2)];
    for(int i = 0; i<k; ++i){
        // https://www.gnu.org/software/gsl/doc/html/linalg.html
        // Examples
    }


//    for(int point = 2; point < k; ++point ){
//        for (int degree = 2; degree < k; ++degree){
//            printf("%10.6f ", coeff((double)point, degree));
//        }
//        printf("= %10.6f\n", f(point));
//    }

    return coeffs;

}

double final_u(const double x, const double *coeffs){
    assert((sizeof(coeffs)/sizeof(coeffs[0])) == k); // should be exactly k coefficients
    double sum = 0;
    for(int i = 0; i < k; ++i){
        sum = sum + coeffs[i]*peano_polynomial(i+1, x); // coeffs[i] stays for alpha(i+1)
    }
    return sum;
}

int main() {
    double* coeffs = solve_linear_system();
    for(int point = a; point <= b; ++point){
        printf("u(%d) = %15.10f\n", point, final_u((double) point, coeffs));
    }
    free(coeffs);
    return 0;
}