/*
 * Project:  This file defines various functions to fit azimuthally averaged data
 *           to a Fermi density distribution.
 *
 * Author:   Pedro M Duarte 2011-01
 * 
 */

#include "fits/fits.h"

extern bool VERBOSE;
using namespace std;

/* Model to evaluate the Fermi azimuthal average.
 * This is just the 2D Fermi with a radiux insted of (x,y)
 */
double
fermiAzimuth_model (double dist, const gsl_vector * v)
{
  double n0 = gsl_vector_get (v, 0);
  double BetaMu = gsl_vector_get (v, 1);
  double r = gsl_vector_get (v, 2);
  double B = gsl_vector_get (v, 3);
  double mx = gsl_vector_get (v, 4);

  return
    n0 / f1 (BetaMu) * f1 (BetaMu - fq (BetaMu) * (pow (dist / r, 2))) +
    B + mx * dist;
}

/* Azimuthal-type data after model evaluation
 */
gsl_vector *
fermiAzimuth_eval (gsl_vector ** d, const double fermi_azimuth_fit[5])
{
  gsl_vector *dat = gsl_vector_alloc (d[0]->size);
  int nparams = 5;
  gsl_vector *v = gsl_vector_alloc (nparams);
  for (int e = 0; e < nparams; e++)
    gsl_vector_set (v, e, fermi_azimuth_fit[e]);
  for (unsigned int i = 0; i < d[0]->size; i++)
    {
      double dist = gsl_vector_get (d[0], i);
      gsl_vector_set (dat, i, fermiAzimuth_model (dist, v));
    }
  return dat;
}

/* Model to evaluate the Fermi azimuthal average at T=0.
 * This is just the 2D Fermi with a radiux insted of (x,y)
 */
double
fermiAzimuthZeroT_model (double dist, const gsl_vector * v)
{
  double n0 = gsl_vector_get (v, 0);
  double r = gsl_vector_get (v, 1);
  double B = gsl_vector_get (v, 2);
  double mx = gsl_vector_get (v, 3);
  return
    n0 * pow (std::max (1. - pow (dist / r, 2.), 0.), 2.) + B + mx * dist;
}

/* T=0 Azimuthal-type data after model evaluation
 */
gsl_vector *
fermiAzimuthZeroT_eval (gsl_vector ** d,
			const double fermi_azimuth_fit_zero[4])
{
  gsl_vector *dat = gsl_vector_alloc (d[0]->size);
  int nparams = 4;
  gsl_vector *v = gsl_vector_alloc (nparams);
  for (int e = 0; e < nparams; e++)
    gsl_vector_set (v, e, fermi_azimuth_fit_zero[e]);
  for (unsigned int i = 0; i < d[0]->size; i++)
    {
      double dist = gsl_vector_get (d[0], i);
      gsl_vector_set (dat, i, fermiAzimuthZeroT_model (dist, v));
    }
  return dat;
}


/**** ERROR FUNCTIONS ****/


double
fermi1d_azimuthal_simplex_f (const gsl_vector * v, void *params)
{
  // params contains an two element array of pointers to gsl_vectors
  // 
  // one gsl_vector contains the distances 
  // and the other contains the azimuthal average at that distance.

  //  benchmark start
  //  double start = omp_get_wtime(); 

  unsigned int s = (((gsl_vector **) params)[0])->size;
  double sumsq = 0.;

  gsl_vector *d = ((gsl_vector **) params)[0];
  gsl_vector *az = ((gsl_vector **) params)[1];

  //i is radial === y
  //j is axial  === g

  for (unsigned int i = 0; i < s; i++)
    {
      //printf("Inisde simplex_f loop\n"); 
      double dist = gsl_vector_get (d, i);
      double dat = gsl_vector_get (az, i);
      sumsq += pow (fermiAzimuth_model (dist, v) - dat, 2);
    }

//  end benchmark
//  double end = omp_get_wtime();
//  printf( "time elapsed in fermi errfunc %.5f\n", end-start); 

  return sumsq;
}

double
fermi1d_azimuthal_zero_simplex_f (const gsl_vector * v, void *params)
{
  //benchmark start
//  double start = omp_get_wtime(); 

  unsigned int s = (((gsl_vector **) params)[0])->size;



  double sumsq = 0.;

  gsl_vector *d = ((gsl_vector **) params)[0];
  gsl_vector *az = ((gsl_vector **) params)[1];

  //i is radial === y
  //j is axial  === g

  for (unsigned int i = 0; i < s; i++)
    {
      //printf("Inisde simplex_f loop\n"); 
      double dist = gsl_vector_get (d, i);
      double dat = gsl_vector_get (az, i);
      sumsq += pow (fermiAzimuthZeroT_model (dist, v) - dat, 2);
    }

//  end benchmark
//  double end = omp_get_wtime();
//  printf( "time elapsed in fermi errfunc %.5f\n", end-start); 

  return sumsq;
}


/**** FITTING LOOPS ****/


void
fit1dfermi_azimuthal_neldermead (gsl_vector ** a, double *fit)
{

  double n0 = fit[0];
  double BetaMu = fit[1];
  double r = fit[2];
  double B = fit[3];
  double mx = fit[4];

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function f;

/*Starting Point */
  x = gsl_vector_alloc (5);
  gsl_vector_set (x, 0, n0);
  gsl_vector_set (x, 1, BetaMu);
  gsl_vector_set (x, 2, r);
  gsl_vector_set (x, 3, B);
  gsl_vector_set (x, 4, mx);


/*Set initial step sizes to 1*/
  ss = gsl_vector_alloc (5);
  gsl_vector_set_all (ss, 1.0);

  f.n = 5;
  f.f = &fermi1d_azimuthal_simplex_f;
  f.params = a;


  s = gsl_multimin_fminimizer_alloc (T, 5);

  gsl_multimin_fminimizer_set (s, &f, x, ss);


  size_t iter = 0;
  int status;
  double size;


  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);
      if (status)
	break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-6);

      if (status == GSL_SUCCESS && VERBOSE)
	{
	  printf (" converged to minimum at \n");
	}

      if ((VERBOSE) && iter % 20 == 0)
	{
	  printf
	    ("%5d %10.3e %10.3e %10.3e %10.3e %10.3e f() = %7.3e size = %.5f\n",
	     iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	     gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3),
	     gsl_vector_get (s->x, 4), s->fval, size);
	}
    }
  while (status == GSL_CONTINUE && iter < 1000);

  gsl_vector_free (x);
  gsl_vector_free (ss);
  gsl_multimin_fminimizer_free (s);

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);
  fit[4] = FIT (4);
  return;
}


void
fit1dfermi_azimuthal_zero_neldermead (gsl_vector ** a, double *fit)
{

  double n0 = fit[0];
  double r = fit[1];
  double B = fit[2];
  double mx = fit[3];

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function f;

/*Starting Point */
  x = gsl_vector_alloc (4);
  gsl_vector_set (x, 0, n0);
  gsl_vector_set (x, 1, r);
  gsl_vector_set (x, 2, B);
  gsl_vector_set (x, 3, mx);


/*Set initial step sizes to 1*/
  ss = gsl_vector_alloc (4);
  gsl_vector_set_all (ss, 1.0);

  f.n = 4;
  f.f = &fermi1d_azimuthal_zero_simplex_f;
  f.params = a;


  s = gsl_multimin_fminimizer_alloc (T, 4);

  gsl_multimin_fminimizer_set (s, &f, x, ss);


  size_t iter = 0;
  int status;
  double size;


  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);
      if (status)
	break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-6);

      if (status == GSL_SUCCESS && VERBOSE)
	{
	  printf (" converged to minimum at \n");
	}

      if ((VERBOSE) && iter % 20 == 0)
	{
	  printf
	    ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3e size = %.5f\n",
	     iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	     gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3),
	     s->fval, size);
	}
    }
  while (status == GSL_CONTINUE && iter < 1000);

  gsl_vector_free (x);
  gsl_vector_free (ss);
  gsl_multimin_fminimizer_free (s);

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);
  return;
}
