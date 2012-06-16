/*
 * Project:  Implementations of fourier transforms for image filtering
 *
 * Author:   Pedro M Duarte 2011-03
 * 
 */

#include "fourier.h"
#include "utils/utils.h"

extern bool VERBOSE;

using namespace std;


gsl_matrix *
fft2d (gsl_matrix * m)
{

  unsigned int s1 = m->size1;
  unsigned int s2 = m->size2;

  gsl_matrix *FT = gsl_matrix_alloc (s1, s2);

  gsl_complex F[s1][s2];
  gsl_complex c;

  for (unsigned int u = 0; u < s1; u++)
    {
      for (unsigned int v = 0; v < s2; v++)
	{
	  F[u][v] = gsl_complex_polar (0., 0.);
	  for (unsigned int x = 0; x < s1; x++)
	    {
	      for (unsigned int y = 0; y < s2; y++)
		{

		  c =
		    gsl_complex_polar (gsl_matrix_get (m, x, y),
				       -2 * M_PI * (u * x / s1 + v * y / s2));
		  F[u][v] = gsl_complex_add (F[u][v], c);
	    }}
	  gsl_matrix_set (FT, u, v, gsl_complex_abs (F[u][v]));
	}
      cout << u << "\t";
    }

  return FT;

}
