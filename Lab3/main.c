#include <stdio.h>
#include <math.h>
#include <assert.h>

#define a 0
#define b 19
#define k 19
#define h 1 // (b-a)/h
#define x_a 0
#define x_b 0


double p(double point){
    return point*point + 1;
}

double q(double point){
    return point*point + cos(point);
}

double f(double point){
    return point*point*sin(point);
}

double peano_polynomial(int degree, double point){
    assert(point >= a && point <= b);
    assert(degree > 0);
    switch (degree) {
        case 1:  return (b - point)/(b - a);
        case 2:  return (point - a)/(b - a);
        case 3:  return ((b - point)*(point - a))/((b-a)*(b-a));
        default: return ((b - point)*(point - a)*pow((2*point-a-b), degree-3))/pow((b-a), degree-1);
    }
}

double first_derivative(double (*f)(int, double), int degree, double point){
    double next_point = point + h;
    double previous_point = point - h;
    return (f(degree, next_point) - f(degree, previous_point))/(2*h);
}

double second_derivative(double (*f)(int, double), int degree, double point){
    double last_point = point + h;
    double previous_point = point - h;
    return (f(degree, previous_point) + f(degree, last_point) - 2*f(degree, point))/(h*h);
}

double coeff(double (*p)(double), double (*q)(double), double point, int degree){
    return second_derivative(&peano_polynomial, degree, point)           +
           p(point) * first_derivative(&peano_polynomial, degree, point) +
           q(point) * peano_polynomial(degree, point);
}


int main() {
    for(int point = 2; point < k; ++point ){
        for (int degree = 2; degree < k; ++degree){
            printf("%10.6f ", coeff(&p, &q, (double)point, degree));
        }
        printf("= %10.6f\n", f(point));
    }

    return 0;
}