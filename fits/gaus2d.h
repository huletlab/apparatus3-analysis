/*
 * Project:  This file defines various functions to fit data to a 2D Gaussian
 *           Several methods are implemented, see the function prototypes for
 *           a list 
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

//  Fitting methods:
void fit2dgaus_neldermead (gsl_matrix *m, double *fit);
void fit2dgaus( gsl_matrix * m, double *fit);
void fit2dgaus_err( gsl_matrix * m, double *fit, double *err);
void fit2dgaus_no_offset( gsl_matrix * m, double *fit);

// Shorthand syntax for getting fit results and fit errors
#define FIT(i) gsl_vector_get( s->x , i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

/* Model to evaluate the 2D Gaussian
 *
 */
double gaus2d_model( double i, double j, const gsl_vector *v){
  // In this model width is 1/e
  double cx = gsl_vector_get (v, 0);
  double wx = gsl_vector_get (v, 1);
  double cy = gsl_vector_get (v, 2);
  double wy = gsl_vector_get (v, 3);
  double A = gsl_vector_get (v, 4);
  double B; 
  if ( v->size == 6 ) B=gsl_vector_get (v, 5);
  else if ( v->size == 5) B=0.; 
  else { printf(" ERROR --> Number of parameters in gaus2d model is invalid "); exit(1); } 
  
  return    B +
	    A * exp (-1.* (pow ((i - cx) / wx, 2.) + pow ((j - cy) / wy, 2.)));
}

/* Matrix data with 2D Gaussian evaluation
 *
 */
gsl_matrix *gaus2d_eval( const gsl_matrix  *d, const double gaus_fit[6], const bool offset = true){
  gsl_matrix *eval = gsl_matrix_alloc( d->size1, d->size2 ) ; 
  int nparams = offset==true ? 6 :5;
   
  gsl_vector *v = gsl_vector_alloc( nparams );
//  fprintf( stderr, "\nEvaulating 2D Gaussian fit results. nparams=%d\n",nparams);  
  for ( int e=0; e < nparams; e++){
  gsl_vector_set(v, e, gaus_fit[e] ) ; 
//  fprintf( stderr, "gaus2d_fit[%d] = %f\n", e, gaus_fit[e]); 
    }  
  
  for (unsigned int i =0; i< d->size1; i++)
   { 
     for ( unsigned int j=0; j<d->size2; j++)
     {
       gsl_matrix_set( eval, i, j, gaus2d_model( i, j, v )); 
    }
  }
  return eval; 
}


/* Matrix data with residuals
 *
 */
gsl_matrix *gaus2d_residual( const gsl_matrix  *d, const double gaus_fit[6], const bool offset = true){
  gsl_matrix *residual = gsl_matrix_alloc( d->size1, d->size2 ) ; 
  int nparams = offset==true ? 6 :5;
   
  gsl_vector *v = gsl_vector_alloc( nparams );
//  fprintf( stderr, "\nEvaulating 2D Gaussian fit results. nparams=%d\n",nparams);  
  for ( int e=0; e < nparams; e++){
  gsl_vector_set(v, e, gaus_fit[e] ) ; 
//  fprintf( stderr, "gaus2d_fit[%d] = %f\n", e, gaus_fit[e]); 
    }  
  
  for (unsigned int i =0; i< d->size1; i++)
   { 
     for ( unsigned int j=0; j<d->size2; j++)
     {
       gsl_matrix_set( residual, i, j,  gsl_matrix_get( d, i, j) - gaus2d_model( i, j, v )); 
    }
  }
  return residual; 
}


void make_gaus2d_inspect( gsl_matrix *c, const double gaus2d_fit[6], const char *prefix){
  string datfile (prefix);
  datfile += "_gaus2ddat.ascii"; 
  string fitfile (prefix); 
  fitfile += "_gaus2dfit.ascii";
 
  save_gsl_matrix_ASCII ( c , datfile);  
  gsl_matrix *fit2d = gaus2d_eval (c, gaus2d_fit);
  save_gsl_matrix_ASCII (fit2d, fitfile);
 
  stringstream inspectstr;
  inspectstr << "inspect2d_ascii.py ";
  inspectstr << datfile;
  inspectstr << " ";
  inspectstr << fitfile;
  inspectstr << " ";
  inspectstr << floor (gaus2d_fit[2]);
  inspectstr << " ";
  inspectstr << floor (gaus2d_fit[0]);
  inspectstr << " ";
  inspectstr << prefix;
  inspectstr << "_gaus";
  //cerr << endl << inspectstr.str () << endl;
  system (inspectstr.str ().c_str ());

  remove ( datfile.c_str());
  remove ( fitfile.c_str());
  return; 
} 
 

/* Error function used in the implementation of the Nelder-Mead
 * fitting algorithm.
 *
 */ 
double
gaus2d_simplex_f (const gsl_vector * v, void *params)
{
  unsigned int s1 = ((gsl_matrix *) params)->size1;
  unsigned int s2 = ((gsl_matrix *) params)->size2;

  double sumsq = 0.;

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double model = gaus2d_model( (double) i, (double) j, v); 
	  double dat = gsl_matrix_get ((gsl_matrix *) params, i, j);
	  sumsq += pow (model - dat, 2);
	}
    }

  return sumsq;
}

/* Error function used in the implementation of the Nelder-Mead
 * fitting algorithm with no offset
 * 
 */ 
double
gaus2d_no_offset_simplex_f (const gsl_vector * v, void *params)
{
  unsigned int s1 = ((gsl_matrix *) params)->size1;
  unsigned int s2 = ((gsl_matrix *) params)->size2;

  double sumsq = 0.;

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double model = gaus2d_model( (double) i, (double) j, v); 
	  double dat = gsl_matrix_get ((gsl_matrix *) params, i, j);
	  sumsq += pow (model - dat, 2);
	}
    }

  return sumsq;
}

/* Nelder-Mead algorithm with no offset
 *
 */
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

/* Nelder-Mead algorithm
 * 
 */ 
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


/* Error vector for the implementation of the Levenberg-Marquardt fitting
 * algorithm
 * See fit2dgaus
 *
 */ 

int
gaus2d_f (const gsl_vector * x, void *data, gsl_vector * f)
{
  unsigned int s1 = ((gsl_matrix *) data)->size1;
  unsigned int s2 = ((gsl_matrix *) data)->size2;

  size_t ii = 0;

  //#pragma omp parallel for
  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double model = gaus2d_model( (double) i, (double) j, x); 
	  double dat = gsl_matrix_get ((gsl_matrix *) data, i, j);
	  gsl_vector_set (f, ii, (model - dat));
	  ii++;
	}
    }

  return GSL_SUCCESS;
}


/* Jacobian matrix for the implementation of the Levenerg-Marquardt fitting
 * algorithm
 */
int
gaus2d_df (const gsl_vector * x, void *data, gsl_matrix * J)
{

  unsigned int s1 = ((gsl_matrix *) data)->size1;
  unsigned int s2 = ((gsl_matrix *) data)->size2;

  double cx = gsl_vector_get (x, 0);
  double wx = gsl_vector_get (x, 1);
  double cy = gsl_vector_get (x, 2);
  double wy = gsl_vector_get (x, 3);
  double A = gsl_vector_get (x, 4);
//  double B = gsl_vector_get (x, 5);

  size_t ii = 0;

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double E =
	    exp (-1 * (pow ((i - cx) / wx, 2) + pow ((j - cy) / wy, 2)));
	  double df_dcx = A * E * 2. * (i - cx) / pow (wx, 2);
	  double df_dwx = A * E * pow (i - cx, 2) * 2.0 / pow (wx, 3);
	  double df_dcy = A * E * 2. * (j - cy) / pow (wy, 2);
	  double df_dwy = A * E * pow (j - cy, 2) * 2.0 / pow (wy, 3);
	  double df_dA = E;
	  double df_dB = 1;

	  gsl_matrix_set (J, ii, 0, df_dcx);
	  gsl_matrix_set (J, ii, 1, df_dwx);
	  gsl_matrix_set (J, ii, 2, df_dcy);
	  gsl_matrix_set (J, ii, 3, df_dwy);
	  gsl_matrix_set (J, ii, 4, df_dA);
	  gsl_matrix_set (J, ii, 5, df_dB);

	  ii++;
	}
    }
  return GSL_SUCCESS;
}

/* Combines the error vector and the Jacobian matrix.  Used in the 
 * Levenberg-Marquardt algorithm. 
 * See fit2dgaus
 *
 */ 
int
gaus2d_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  gaus2d_f (x, data, f);
  gaus2d_df (x, data, J);

  return GSL_SUCCESS;
}



/* Function to print the status of the Levenberg-Marquardt solver
 *
 */
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



/* Levenberg-Marquardt fitting algorithm
 * 
 */ 
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


/* Levenberg-Marquardt fitting algorithm
 * This one also returns the fit errors
 * 
 */ 
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
