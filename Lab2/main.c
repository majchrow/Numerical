/*
 * @author Dawid Majchrowski
 * @date   12.11.2019
 * Plot curve f(x) = Int_{a}^{x} f(x) dx | x in interval [a,b]
*/

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <assert.h>
#include "dislin.h"

// Constant definitions
#define a    0.                       // lower bound
#define b    50.                      // upper bound
#define n    1000                     // grid size for plotting

// Macros
#define MAX(A,B) (A) > (B) ? (A) : (B)
#define MIN(A,B) (A) > (B) ? (B) : (A)

// Method definitions
void plot_result(int N);
double f(double x, void *params);
double int_f(double lower_bound, double upper_bound);


// Method implementations
double f(double x, void *params) {
    assert(a <= x && x <= b);
    double alpha = *(double *) params;
    double f = alpha * exp(sin(x)); // alpha = 1
    return f;
}

double int_f(const double lower_bound, const double upper_bound){ // int_{lower_bound}^{upper_bound} f(x) dx
    assert(lower_bound <= upper_bound);
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
    double result, error;
    double epsabs = 0;
    double epsrel = 1e-7;
    size_t limit = 1000;
    gsl_function F;
    F.function = &f;
    double aplha = 1.0;
    F.params = &aplha;
    gsl_integration_qags (&F, lower_bound, upper_bound, epsabs, epsrel, limit, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);
    return result;
}

void plot_result(const int N){
    float *arguments = calloc(N, sizeof(double));
    float *values    = calloc(N, sizeof(double));
    double step = (b-a)/(N-1);
    int ic =   intrgb (0.95,0.95,0.95);
    double x_l = a - 5.*step;             // x_left
    double x_r = b + 5.*step;             // x_right
    double y_b = 0.; // y_bottom
    double y_t = 0.; // y_top

    for (int i = 0; i < N; ++i){
        arguments[i] = a + i*step;
        values[i] = int_f(a, arguments[i]);
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
    char title_desc[4096];
    sprintf(title_desc, "f(t) = Int_{%1.1f}^{t}(e^(sin(x)))dx, t e [%1.1f, %1.1f]", a, a, b);
    titlin (title_desc, 1);
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
           "number of grid points = %d\n", a, b, n);
    plot_result(n+1);
    return 0;
}
