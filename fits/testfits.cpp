/*
 * Project:  Testing fitting routines
 *
 * Author:   Pedro M Duarte 2011-02
 * 
 */

#include "fits.h"
#include "utils/utils.h"
#include "funcs/funcs.h"

using namespace std;

bool VERBOSE;

int
main (int argc, char **argv)
{

  VERBOSE = false;

  string file ("0649.fluor");
  gsl_matrix *img = ReadFluorImg_gsl_matrix (file);

/*  gsl_vector *row = gsl_vector_alloc (img->size2);
  for (unsigned int j = 0; j < img->size2; j++)
    {
      gsl_vector_set (row, j, gsl_matrix_get (img, 374, j));
      //    cout << gsl_vector_get (row, j) << "\t" << endl;
    }

  double gaus_fit[4] = { 481.15, 23.82, 51.5, .1 };
  fit1dgaus_neldermead (row, gaus_fit);
  fit1dgaus (row, gaus_fit); */

  gsl_matrix *c = autocropImage (img, 6.);
  string fileout2 ("0649c.fluor");
  save_gsl_matrix_ASCII (c, fileout2);

  unsigned int wi, wj, ci, cj;
  double max;
  findpeak (c, &ci, &cj, &max);
  findFWHM (c, &wi, &wj);

  gsl_vector *max_row = gsl_vector_alloc (c->size2);
  for (unsigned int j = 0; j < c->size2; j++)
    {
      gsl_vector_set (max_row, j, gsl_matrix_get (c, ci, j));
    }

  gsl_vector *max_col = gsl_vector_alloc (c->size1);
  for (unsigned int i = 0; i < c->size1; i++)
    {
      gsl_vector_set (max_col, i, gsl_matrix_get (c, i, cj));
    }
  double gaus_fit_i[4] = { (double) cj, (double) wj, max, 0.1 };
  fit1dgaus (max_row, gaus_fit_i);

  double gaus_fit_j[4] = { (double) ci, (double) wi, max, 0.1 };
  fit1dgaus (max_col, gaus_fit_j);

  double gaus2d_fit[6] =
    { gaus_fit_i[1], gaus_fit_i[2], gaus_fit_j[1], gaus_fit_j[2], max, 0.1 };

//  double gaus2d_fit[6] = { 62., 20, 49., 17.4, 51.5, .1 };
  fit2dgaus_neldermead (c, gaus2d_fit);
  fit2dgaus (c, gaus2d_fit);



  cout << "testing ..." << endl;

  return 0;
}
