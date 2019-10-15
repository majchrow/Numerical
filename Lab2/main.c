#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "dislin.h"

#define ABS(X) (X)>0?(X):-(X)


double f (double x, void *params) {
  double f = exp(sin(x));
  return f;
}

int main ()
{
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

  double lower_bound = 0;
  double upper_bound = 0;
  double result, error;
  double epsabs = 0;
  double epsrel = 1e-7;
  size_t limit = 1000;
  gsl_function F;
  F.function = &f;

  for(upper_bound = 0; ABS(upper_bound - 9) > epsrel; upper_bound+=0.2){
  	gsl_integration_qags (&F, lower_bound, upper_bound, epsabs, epsrel, limit, workspace, &result, &error);  	
  	printf ("Upper bound=%.2f, result = %.18f\n", upper_bound, result);
  }	

  gsl_integration_workspace_free(workspace);
  return 0;
}
