/*
 * Project:  Model functions for 2D gaussian fit
 *
 * Author:   Pedro M Duarte 2011-01
 * 
 */


#include<math.h>

double
gaus2d_simplex_f (const gsl_vector * v, void *params)
{
  unsigned int s1 = ((gsl_matrix *) params)->size1;
  unsigned int s2 = ((gsl_matrix *) params)->size2;

  double cx = gsl_vector_get (v, 0);
  double wx = gsl_vector_get (v, 1);
  double cy = gsl_vector_get (v, 2);
  double wy = gsl_vector_get (v, 3);
  double A = gsl_vector_get (v, 4);
  double B = gsl_vector_get (v, 5);

  double sumsq = 0.;

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double model =
	    B +
	    A * exp (-1.* (pow ((i - cx) / wx, 2.) + pow ((j - cy) / wy, 2.)));
	  double dat = gsl_matrix_get ((gsl_matrix *) params, i, j);
	  sumsq += pow (model - dat, 2);
	}
    }

  return sumsq;
}

double
gaus2d_no_offset_simplex_f (const gsl_vector * v, void *params)
{
  unsigned int s1 = ((gsl_matrix *) params)->size1;
  unsigned int s2 = ((gsl_matrix *) params)->size2;

  double cx = gsl_vector_get (v, 0);
  double wx = gsl_vector_get (v, 1);
  double cy = gsl_vector_get (v, 2);
  double wy = gsl_vector_get (v, 3);
  double A = gsl_vector_get (v, 4);

  double sumsq = 0.;

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double model =
	    A * exp (-1.* (pow ((i - cx) / wx, 2.) + pow ((j - cy) / wy, 2.)));
	  double dat = gsl_matrix_get ((gsl_matrix *) params, i, j);
	  sumsq += pow (model - dat, 2);
	}
    }

  return sumsq;
}



int
gaus2d_f (const gsl_vector * x, void *data, gsl_vector * f)
{
  unsigned int s1 = ((gsl_matrix *) data)->size1;
  unsigned int s2 = ((gsl_matrix *) data)->size2;

  double cx = gsl_vector_get (x, 0);
  double wx = gsl_vector_get (x, 1);
  double cy = gsl_vector_get (x, 2);
  double wy = gsl_vector_get (x, 3);
  double A = gsl_vector_get (x, 4);
  double B = gsl_vector_get (x, 5);

  size_t ii = 0;

  //#pragma omp parallel for
  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double model =
	    B +
	    A * exp (-1. * (pow ((i - cx) / wx, 2.) + pow ((j - cy) / wy, 2.)));
	  double dat = gsl_matrix_get ((gsl_matrix *) data, i, j);
	  gsl_vector_set (f, ii, (model - dat));
	  ii++;
	}
    }

  return GSL_SUCCESS;
}

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

int
gaus2d_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  gaus2d_f (x, data, f);
  gaus2d_df (x, data, J);

  return GSL_SUCCESS;
}
