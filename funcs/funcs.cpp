/*
 * Project:  Wrappers for GSL Fermi-Dirac integrals
 *
 * Author:   Pedro M Duarte 2011-01
 * 
 */

#include "funcs.h"
#include "math.h"

extern bool VERBOSE;

using namespace std;


double
f1 (double x)
{
  return gsl_sf_fermi_dirac_1 (x);
}


double
f2 (double x)
{
  return gsl_sf_fermi_dirac_2 (x);
}

double
f32 (double x)
{
  return gsl_sf_fermi_dirac_3half (x);
}

double
fm1 (double x)
{
  return gsl_sf_fermi_dirac_m1 (x);
}

double
f0 (double x)
{
  return gsl_sf_fermi_dirac_0 (x);
}

double
fq (double x)
{
  return (1 + exp (x)) / exp (x) * log (1 + exp (x));
}
