/*
 * Project:  Model functions for 2D gaussian fit
 *
 * Author:   Pedro M Duarte 2011-01
 * 
 */


#include<math.h>

double
gaus1d_simplex_f (const gsl_vector * v, void *params)
{
  unsigned int s1 = ((gsl_vector *) params)->size;

  double cx = gsl_vector_get (v, 0);
  double wx = gsl_vector_get (v, 1);
  double A = gsl_vector_get (v, 2);
  double B = gsl_vector_get (v, 3);

  double sumsq = 0.;

  for (unsigned int i = 0; i < s1; i++)
    {
	  double model = B + A * exp (-1.0 * (pow ((i - cx) / wx, 2) ));
	  double dat = gsl_vector_get ((gsl_vector *) params, i);
	  sumsq += pow (model - dat, 2);
    }

  return sumsq;
}



int
gaus1d_f (const gsl_vector * x, void *data, gsl_vector * f)
{
  unsigned int s1 = ((gsl_vector *) data)->size;

  double cx = gsl_vector_get (x, 0);
  double wx = gsl_vector_get (x, 1);
  double A = gsl_vector_get (x, 2);
  double B = gsl_vector_get (x, 3);

  size_t ii = 0;

  for (unsigned int i = 0; i < s1; i++)
    {
	  double model = B + A * exp (-1.0 * (pow ((i - cx) / wx, 2) ));
	  double dat = gsl_vector_get ((gsl_vector *) data, i);
	  gsl_vector_set (f, ii, model - dat);
	  ii++;
    }

  return GSL_SUCCESS;
}

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

int
gaus1d_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  gaus1d_f (x, data, f);
  gaus1d_df (x, data, J);

  return GSL_SUCCESS;
}
