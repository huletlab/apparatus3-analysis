/*
 * Project:  This file defines various functions to fit data to a 1D Gaussian
 *           Several methods are implemented,  see the function prototypes for 
 *           a list.  
 * 
 *           
 *
 * Author:   Pedro M Duarte 2011-01
 * 
 */


#include <iostream>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>

extern bool VERBOSE;
using namespace std;


// Fitting methods
void fit1dgaus_neldermead (gsl_vector *m, double *fit);
void fit1dgaus( gsl_vector * m, double *fit);

// Shorthand syntax for getting fit results and fit errors. 
#define FIT(i) gsl_vector_get( s->x , i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))


/* Model to evaluate the 1D Gaussian
 *
 */
double gaus1d_model( double x, const  gsl_vector *v){
  // In this model width is 1/e 
  // The model parameters are 
//  double cx = gsl_vector_get (v, 0);
//  double wx = gsl_vector_get (v, 1);
//  double A = gsl_vector_get (v, 2);
//  double B = gsl_vector_get (v, 3);

//  return B + A * exp ( -1.0 * ( pow (( x - cx) / wx, 2) ));
  return gsl_vector_get(v,3) + gsl_vector_get(v,2) * exp ( -1.0 * ( pow (( x - gsl_vector_get(v,0) ) / gsl_vector_get(v,1), 2) ) );  
}


/* Vector data with 1D Gaussian evaluation
 * 
 * v is the data vector, it is only provided so that the output
 * vector has the same length
 *
 */
gsl_vector *gaus1d_eval(const gsl_vector *d, const double gaus_fit[4]){
  gsl_vector *eval = gsl_vector_alloc( d->size); 
  int nparams = 4;
  gsl_vector *v = gsl_vector_alloc( nparams); 
  for ( int e=0; e<nparams; e++){
    gsl_vector_set(v, e, gaus_fit[e]);
  } 
  for ( unsigned int i =0;  i < v->size ; i++){
    gsl_vector_set( eval, i, gaus1d_model( i, v));
  }
  return eval;
}

   

/* Error function used in the implementation of the Nelder-Mead
   fitting algorithm.
   See fit1dgaus_neldermead   
 *
 */ 
double
gaus1d_simplex_f (const gsl_vector * v, void *params)
{
  unsigned int s1 = ((gsl_vector *) params)->size;
  double sumsq = 0.;

  for (unsigned int i = 0; i < s1; i++)
    {
	  double model = gaus1d_model( i, v) ;
	  double dat = gsl_vector_get ((gsl_vector *) params, i);
	  sumsq += pow (model - dat, 2);
    }

  return sumsq;
}


/* Nelder-Mead fitting algorithm
 *
 */ 
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



/* Error vector for the implementation of the Levenberg-Marquardt fitting
   algorithm
   See fit1dgaus
 * 
 */
int
gaus1d_f (const gsl_vector * x, void *data, gsl_vector * f)
{
  unsigned int s1 = ((gsl_vector *) data)->size;
  size_t ii = 0;

  for (unsigned int i = 0; i < s1; i++)
    {
	  double model = gaus1d_model( i, x) ;
	  double dat = gsl_vector_get ((gsl_vector *) data, i);
	  gsl_vector_set (f, ii, model - dat);
	  ii++;
    }

  return GSL_SUCCESS;
}



/* Jacobiam matrix for the implementation of the Levenberg-Marquardt fitting
   algorithm
   See fit1dgaus
 * 
 */
int
gaus1d_df (const gsl_vector * x, void *data, gsl_matrix * J)
{

  unsigned int s1 = ((gsl_vector *) data)->size;

  double cx = gsl_vector_get (x, 0);
  double wx = gsl_vector_get (x, 1);
  double A = gsl_vector_get (x, 2);
//  double B = gsl_vector_get (x, 5);

  size_t ii = 0;

  for (unsigned int i = 0; i < s1; i++)
    {
	  double E =  exp (-1.0 * (pow ((i - cx) / wx, 2) ));
	  double df_dcx = A * E * 2. * (i - cx) / pow (wx, 2);
	  double df_dwx = A * E * pow (i - cx, 2) * 2.0 / pow (wx, 3);
	  double df_dA = E;
	  double df_dB = 1;

	  gsl_matrix_set (J, ii, 0, df_dcx);
	  gsl_matrix_set (J, ii, 1, df_dwx);
	  gsl_matrix_set (J, ii, 2, df_dA);
	  gsl_matrix_set (J, ii, 3, df_dB);

	  ii++;
    }
  return GSL_SUCCESS;
}


/* Combines the error vector and the Jacobian matrix.  Used in the 
 * Levenberg-Marquardt algorithm. 
 * See fit1dgaus
 *
 */ 
int
gaus1d_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  gaus1d_f (x, data, f);
  gaus1d_df (x, data, J);

  return GSL_SUCCESS;
}


/* Function to print the status of the Levenberg-Marquardt solver
 *
 */
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


/* Levenberg-Marquardt fitting algorithm
 * 
 */ 
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
