#ifndef FITS_HEADER_FILE
#define FITS_HEADER_FILE

#include <omp.h>

#include <iostream>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>

#include "funcs/funcs.h"
#include "utils/utils.h"


// Shorthand syntax for getting fit results and fit errors. 
#define FIT(i) gsl_vector_get( s->x , i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))


/* 
 * Functions in gaus1d.cpp that will be exported:
 *
 */
double gaus1d_model (double x, const gsl_vector * v);
gsl_vector *gaus1d_eval (const gsl_vector * d, const double gaus_fit[4]);

void fit1dgaus_neldermead (gsl_vector * m, double *fit);
void fit1dgaus (gsl_vector * m, double *fit);


/* 
 * Functions in gaus2d.cpp that will be exported:
 *
 */
double gaus2d_model (double i, double j, const gsl_vector * v);
double mottGaus2d_model (double i, double j, const gsl_vector * v);
gsl_matrix *gaus2d_eval (const gsl_matrix * d, const double gaus_fit[6],
			 const bool offset = true);
void gaus2d_eval_Azimuth (const double gaus2dfit[6], string prefix);
gsl_matrix *gaus2d_residual (const gsl_matrix * d, const double gaus_fit[6],
			     const bool offset = true);
void make_gaus2d_inspect (gsl_matrix * c, const double gaus2d_fit[6],
			  const char *prefix);

void fit2dgaus_neldermead (gsl_matrix * m, double *fit);
void fit2dmottgaus_neldermead (gsl_matrix * m, double *fit);
void fit2dgaus (gsl_matrix * m, double *fit);
void fit2dgaus_err (gsl_matrix * m, double *fit, double *err);
void fit2dgaus_no_offset (gsl_matrix * m, double *fit);

/* This one takes care of doing all the initial guessing */
void Fit2DGaus_High_Level (gsl_matrix * m, double *fit, double *fit_err,
			   string prefix);
void Gaus2DGuess (gsl_matrix * m, double *guess, string prefix,
		  bool save_matrices);


/* 
 * Functions in fermi2d.cpp that will be exported:
 *
 */
double fermi2d_model (double i, double j, const gsl_vector * v);
gsl_matrix *fermi2d_eval (const gsl_matrix * d, const double fermi_fit[7]);
void fermi2d_eval_Azimuth (const double fermi_fit[7], string prefix);
void make_fermi2d_inspect (gsl_matrix * c, const double fermi2d_fit[6],
			   const char *prefix);
void make_fermi2d_gaus2d_inspect (gsl_matrix * c, const double fermi2d_fit[7],
				  const double gaus2d_fit[6],
				  const char *prefix);


void fit2dfermi_neldermead (gsl_matrix * m, double *fit);
void fit1dfermi_neldermead (gsl_vector * m, double *fit);

/*
 * Functions in fermiAzimuth.cpp that will be exported:
 *
 */
double fermiAzimuth_model (double dist, const gsl_vector * v);
gsl_vector *fermiAzimuth_eval (gsl_vector * d,
			       const double fermi_azimuth_fit[5]);
double fermiAzimuthZeroT_model (double dist, const gsl_vector * v);
gsl_vector *fermiAzimuthZeroT_eval (gsl_vector * d,
				    const double fermi_azimuth_fit_zero[4]);
void fit1dfermi_azimuthal_neldermead (gsl_vector ** a, double *fit);
void fit1dfermi_azimuthal_zero_neldermead (gsl_vector ** a, double *fit);






#endif
