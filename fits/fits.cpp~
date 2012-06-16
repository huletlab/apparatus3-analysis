/*
 * Project:  Implementations of general fitting routines
 *
 * Author:   Pedro M Duarte 2011-01
 * 
 */

#include "fits.h"
#include "funcs/funcs.h"

#include <math.h>
#include "gaus1d.h"
#include "gaus2d.h"
#include "fermi2d.h"

extern bool VERBOSE;

using namespace std;

void
print_state_1d (size_t iter, gsl_multifit_fdfsolver * s)
{
  if (VERBOSE)
    {
      printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f % 15.8f  "
	      "|f(x)| = %g\n",
	      iter,
	      gsl_vector_get (s->x, 0),
	      gsl_vector_get (s->x, 1),
	      gsl_vector_get (s->x, 2),
	      gsl_vector_get (s->x, 3), gsl_blas_dnrm2 (s->f));
    }

}

#define FIT(i) gsl_vector_get( s->x , i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

void
fit1dgaus_neldermead (gsl_vector * m, double *fit)
{

  double cx = fit[0];
  double wx = fit[1];
  double A = fit[2];
  double B = fit[3];

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function f;

/*Starting Point */
  x = gsl_vector_alloc (4);
  gsl_vector_set (x, 0, cx);
  gsl_vector_set (x, 1, wx);
  gsl_vector_set (x, 2, A);
  gsl_vector_set (x, 3, B);


/*Set initial step sizes to 1*/
  ss = gsl_vector_alloc (4);
  gsl_vector_set_all (ss, 1.0);

  f.n = 4;
  f.f = &gaus1d_simplex_f;
  f.params = m;

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
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS && VERBOSE)
	{
	  printf (" converged to minimum at \n");
	}

      if (VERBOSE)
	{
	  printf
	    ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3e size = %.3f\n",
	     iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	     gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3), s->fval,
	     size);
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


void
fit1dgaus (gsl_vector * m, double *fit)
{

  double cx = fit[0];
  double wx = fit[1];
  double A = fit[2];
  double B = fit[3];


  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int iter = 0;
  const size_t n = (m->size);
  const size_t p = 4;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  gsl_multifit_function_fdf f;


  double x_init[4] = { cx, wx, A, B };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);


  f.f = &gaus1d_f;
  f.df = &gaus1d_df;
  f.fdf = &gaus1d_fdf;
  f.n = n;
  f.p = p;
  f.params = m;

  T = gsl_multifit_fdfsolver_lmder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);
  print_state_1d (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
      if (VERBOSE)
	{
	  printf ("status = %s\n", gsl_strerror (status));
	}
      print_state_1d (iter, s);

      if (status)
	break;
      status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);


  {
    double chi = gsl_blas_dnrm2 (s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL (1, chi / sqrt (dof));

    if (VERBOSE)
      {

	printf ("chisq/dof = %g\n", pow (chi, 2.0) / dof);

	printf ("cx     = %.5f +/- %.5f\n", FIT (0), c * ERR (0));
	printf ("wx     = %.5f +/- %.5f\n", FIT (1), c * ERR (1));
	printf ("A      = %.5f +/- %.5f\n", FIT (2), c * ERR (2));
	printf ("B      = %.5f +/- %.5f\n", FIT (3), c * ERR (3));
      }
  }

  if (VERBOSE)
    {
      printf ("status = %s\n", gsl_strerror (status));
    }

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);

  return;
}

void
fit2dgaus_no_offset (gsl_matrix * m, double *fit)
{

  double cx = fit[0];
  double wx = fit[1];
  double cy = fit[2];
  double wy = fit[3];
  double A = fit[4];

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function f;

/*Starting Point */
  x = gsl_vector_alloc (5);
  gsl_vector_set (x, 0, cx);
  gsl_vector_set (x, 1, wx);
  gsl_vector_set (x, 2, cy);
  gsl_vector_set (x, 3, wy);
  gsl_vector_set (x, 4, A);


/*Set initial step sizes to 1*/
  ss = gsl_vector_alloc (5);
  gsl_vector_set_all (ss, 1.0);

  f.n = 5;
  f.f = &gaus2d_no_offset_simplex_f;
  f.params = m;

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
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS && VERBOSE)
	{
	  printf (" converged to minimum at \n");
	}
      if (VERBOSE)
	{
	  printf
	    ("%5d %10.3e %10.3e %10.3e %10.3e %10.3e f() = %7.3e size = %.3f\n",
	     iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	     gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3),
	     gsl_vector_get (s->x, 4), s->fval, size);
	}
    }
  while (status == GSL_CONTINUE && iter < 1000);

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);
  fit[4] = FIT (4);

  gsl_vector_free (x);
  gsl_vector_free (ss);
  gsl_multimin_fminimizer_free (s);

  return;
}

void
fit2dgaus_neldermead (gsl_matrix * m, double *fit)
{

  double cx = fit[0];
  double wx = fit[1];
  double cy = fit[2];
  double wy = fit[3];
  double A = fit[4];
  double B = fit[5];

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function f;

/*Starting Point */
  x = gsl_vector_alloc (6);
  gsl_vector_set (x, 0, cx);
  gsl_vector_set (x, 1, wx);
  gsl_vector_set (x, 2, cy);
  gsl_vector_set (x, 3, wy);
  gsl_vector_set (x, 4, A);
  gsl_vector_set (x, 5, B);


/*Set initial step sizes to 1*/
  ss = gsl_vector_alloc (6);
  gsl_vector_set_all (ss, 1.0);

  f.n = 6;
  f.f = &gaus2d_simplex_f;
  f.params = m;

  s = gsl_multimin_fminimizer_alloc (T, 6);
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
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS && VERBOSE)
	{
	  printf (" converged to minimum at \n");
	}
      if (VERBOSE)
	{
	  printf
	    ("%5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e f() = %7.3e size = %.3f\n",
	     iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	     gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3),
	     gsl_vector_get (s->x, 4), gsl_vector_get (s->x, 5), s->fval,
	     size);
	}
    }
  while (status == GSL_CONTINUE && iter < 1000);

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);
  fit[4] = FIT (4);
  fit[5] = FIT (5);

  gsl_vector_free (x);
  gsl_vector_free (ss);
  gsl_multimin_fminimizer_free (s);

  return;
}


void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  if (VERBOSE)
    {
      printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f % 15.8f % 15.8f %15.8f "
	      "|f(x)| = %g\n",
	      iter,
	      gsl_vector_get (s->x, 0),
	      gsl_vector_get (s->x, 1),
	      gsl_vector_get (s->x, 2),
	      gsl_vector_get (s->x, 3),
	      gsl_vector_get (s->x, 4),
	      gsl_vector_get (s->x, 5), gsl_blas_dnrm2 (s->f));
    }
}


void
fit2dgaus (gsl_matrix * m, double *fit)
{

  double cx = fit[0];
  double wx = fit[1];
  double cy = fit[2];
  double wy = fit[3];
  double A = fit[4];
  double B = fit[5];


  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int iter = 0;
  const size_t n = (m->size1) * (m->size2);
  const size_t p = 6;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  gsl_multifit_function_fdf f;


  double x_init[6] = { cx, wx, cy, wy, A, B };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);


  f.f = &gaus2d_f;
  f.df = &gaus2d_df;
  f.fdf = &gaus2d_fdf;
  f.n = n;
  f.p = p;
  f.params = m;

  T = gsl_multifit_fdfsolver_lmder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);
  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
      if (VERBOSE)
	{
	  printf ("status = %s\n", gsl_strerror (status));
	}
      print_state (iter, s);

      if (status)
	break;
      status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_multifit_covar (s->J, 0.0, covar);


  {
    double chi = gsl_blas_dnrm2 (s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL (1, chi / sqrt (dof));

    if (VERBOSE)
      {
	printf ("chisq/dof = %g\n", pow (chi, 2.0) / dof);

	printf ("cx     = %.5f +/- %.5f\n", FIT (0), c * ERR (0));
	printf ("wx     = %.5f +/- %.5f\n", FIT (1), c * ERR (1));
	printf ("cy     = %.5f +/- %.5f\n", FIT (2), c * ERR (2));
	printf ("wy     = %.5f +/- %.5f\n", FIT (3), c * ERR (3));
	printf ("A      = %.5f +/- %.5f\n", FIT (4), c * ERR (4));
	printf ("B      = %.5f +/- %.5f\n", FIT (5), c * ERR (5));
      }
  }

  if (VERBOSE)
    {
      printf ("status = %s\n", gsl_strerror (status));
    }

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);
  fit[4] = FIT (4);
  fit[5] = FIT (5);


  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);

  return;
}


void
fit2dgaus_err (gsl_matrix * m, double *fit, double *err)
{

  double cx = fit[0];
  double wx = fit[1];
  double cy = fit[2];
  double wy = fit[3];
  double A = fit[4];
  double B = fit[5];


  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int iter = 0;
  const size_t n = (m->size1) * (m->size2);
  const size_t p = 6;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  gsl_multifit_function_fdf f;


  double x_init[6] = { cx, wx, cy, wy, A, B };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);


  f.f = &gaus2d_f;
  f.df = &gaus2d_df;
  f.fdf = &gaus2d_fdf;
  f.n = n;
  f.p = p;
  f.params = m;

  T = gsl_multifit_fdfsolver_lmder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);
  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
      if (VERBOSE)
	{
	  printf ("status = %s\n", gsl_strerror (status));
	}
      print_state (iter, s);

      if (status)
	break;
      status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_multifit_covar (s->J, 0.0, covar);



  double chi = gsl_blas_dnrm2 (s->f);
  double dof = n - p;
  double c = GSL_MAX_DBL (1, chi / sqrt (dof));

  if (VERBOSE)
    {
      printf ("chisq/dof = %g\n", pow (chi, 2.0) / dof);

      printf ("cx     = %.5f +/- %.5f\n", FIT (0), c * ERR (0));
      printf ("wx     = %.5f +/- %.5f\n", FIT (1), c * ERR (1));
      printf ("cy     = %.5f +/- %.5f\n", FIT (2), c * ERR (2));
      printf ("wy     = %.5f +/- %.5f\n", FIT (3), c * ERR (3));
      printf ("A      = %.5f +/- %.5f\n", FIT (4), c * ERR (4));
      printf ("B      = %.5f +/- %.5f\n", FIT (5), c * ERR (5));
    }


  if (VERBOSE)
    {
      printf ("status = %s\n", gsl_strerror (status));
    }

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);
  fit[4] = FIT (4);
  fit[5] = FIT (5);


  err[0] = c * ERR (0);
  err[1] = c * ERR (1);
  err[2] = c * ERR (2);
  err[3] = c * ERR (3);
  err[4] = c * ERR (4);
  err[5] = c * ERR (5);

  /* printf ("\n--- COVARIANCE MATRIX ---\n");
     for (unsigned int i = 0; i < p; i++)
     {
     for (unsigned int j = 0; j < p; j++)
     {
     printf ("%.3e\t\t", gsl_matrix_get (covar, i, j));
     }
     printf ("\n");
     }
     printf("\n"); */

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);

  return;

}

void
fit2dfermi_neldermead (gsl_matrix * m, double *fit)
{

  //i is radial === y
  //j is axial  === g

  double n0 = fit[0];
  double BetaMu = fit[1];
  double ri = fit[2];
  double rj = fit[3];
  double ci = fit[4];
  double cj = fit[5];
  double B = fit[6];

  if (VERBOSE)
    printf (" Nelder-Mead 2D Fermi fit started   ==>  BetaMu= %.3f\n",
	    BetaMu);


  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function f;

/*Starting Point */
  x = gsl_vector_alloc (7);
  gsl_vector_set (x, 0, n0);
  gsl_vector_set (x, 1, BetaMu);
  gsl_vector_set (x, 2, ri);
  gsl_vector_set (x, 3, rj);
  gsl_vector_set (x, 4, ci);
  gsl_vector_set (x, 5, cj);
  gsl_vector_set (x, 6, B);


/*Set initial step sizes to 1*/
  ss = gsl_vector_alloc (7);
  gsl_vector_set_all (ss, 1.0);

  f.n = 7;
  f.f = &fermi2d_simplex_f;
  f.params = m;

  s = gsl_multimin_fminimizer_alloc (T, 7);
  gsl_multimin_fminimizer_set (s, &f, x, ss);


  size_t iter = 0;
  int status;
  double size;

  double tolerance = 1e-2;

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);
      if (status)
	break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tolerance);

      if (status == GSL_SUCCESS && VERBOSE)
	{
	  printf (" converged to minimum at \n");
	}
      if ((VERBOSE) && iter % 20 == 0)
	{
	  printf
	    ("%5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e f() = %7.3e size = %.3f\n",
	     iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	     gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3),
	     gsl_vector_get (s->x, 4), gsl_vector_get (s->x, 5),
	     gsl_vector_get (s->x, 6), s->fval, size);
	}
    }
  while (status == GSL_CONTINUE && iter < 1000);

  fit[0] = FIT (0);
  fit[1] = FIT (1);
  fit[2] = FIT (2);
  fit[3] = FIT (3);
  fit[4] = FIT (4);
  fit[5] = FIT (5);
  fit[6] = FIT (6);

  gsl_vector_free (x);
  gsl_vector_free (ss);
  gsl_multimin_fminimizer_free (s);

  return;
}

void
fit1dfermi_neldermead (gsl_vector * m, double *fit)
{

  double n0 = fit[0];
  double BetaMu = fit[1];
  double r = fit[2];
  double c = fit[3];
  double B = fit[4];

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function f;

/*Starting Point */
  x = gsl_vector_alloc (5);
  gsl_vector_set (x, 0, n0);
  gsl_vector_set (x, 1, BetaMu);
  gsl_vector_set (x, 2, r);
  gsl_vector_set (x, 3, c);
  gsl_vector_set (x, 4, B);


/*Set initial step sizes to 1*/
  ss = gsl_vector_alloc (5);
  gsl_vector_set_all (ss, 1.0);

  f.n = 5;
  f.f = &fermi1d_simplex_f;
  f.params = m;

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
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS && VERBOSE)
	{
	  printf (" converged to minimum at \n");
	}

      if (VERBOSE && iter % 4 == 0)
	{
	  printf
	    ("%5d %10.3e %10.3e %10.3e %10.3e %10.3e f() = %7.3e size = %.3f\n",
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
