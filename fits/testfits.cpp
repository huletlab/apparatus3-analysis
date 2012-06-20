/*
 * Project:  Testing fitting routines
 *
 * Author:   Pedro M Duarte 2011-02
 * 
 */

#include <math.h>
#include "utils/utils.h"
#include "gaus1d.h"
#include "gaus2d.h"
#include "fermi2d.h"
#include "funcs/funcs.h"

using namespace std;

bool VERBOSE;

int
main (int argc, char **argv)
{

  VERBOSE = true;

  // Load data from file 
  string file ("fitexamples/0649.fluor");
  gsl_matrix *img = ReadFluorImg_gsl_matrix (file);

  // Get a row from the data
  gsl_vector *row = gsl_vector_alloc (img->size2);
  for (unsigned int j = 0; j < img->size2; j++)
    {
      gsl_vector_set (row, j, gsl_matrix_get (img, 374, j));
      //    cout << gsl_vector_get (row, j) << "\t" << endl;
    }

  // Do a 1D Gaussian fit on the row
  fprintf (stderr, "\nFitting row #374 with a 1D Gaussian:");
  double gaus_fit[4] = { 481.15, 23.82, 51.5, .1 };

  // Using the Nelder-Mead algorithm
  fprintf (stderr, "\n\n\t...Using Nelder-Mead : \n\t");
  fit1dgaus_neldermead (row, gaus_fit);
  //gsl_vector *nelder1d = gaus1d_eval (row, gaus_fit);
  for (int e = 0; e < 4; e++)
    fprintf (stderr, "[%d]=%f, ", e, gaus_fit[e]);

  // Using the Levenberg-Marquardt algorithm 
  fprintf (stderr, "\n\n\t...Using Levenberg-Marquardt : \n\t");
  fit1dgaus (row, gaus_fit);
  for (int e = 0; e < 4; e++)
    fprintf (stderr, "[%d]=%f, ", e, gaus_fit[e]);


  //OPTIONAL : Crop the matrix for the 2D FIT
  unsigned int roi[4] = { 330, 300, 275, 250 };
  gsl_matrix *c = cropImage_ROI (roi, img);
  //gsl_matrix *c = img;
  string fileout2 ("fitexamples/0649c.fluor");
  save_gsl_matrix_ASCII (c, fileout2);

  // Do a 2D Gaussian fit on the cropped image
  fprintf (stderr, "\n\nFitting cropped matrix with a 2D Gaussian:\n");
  double gaus2d_fit[6] = { 150., 20, 106., 20, 51.5, .1 };


  // Using the Levenberg-Marquardt algorithm
  fprintf (stderr, "\n\n\t...Using Levenberg-Marquardt : \n\t");
  fit2dgaus (c, gaus2d_fit);
  make_gaus2d_inspect (c, gaus2d_fit, "fitexamples/0649lm");

  for (int e = 0; e < 6; e++)
    fprintf (stderr, "[%d]=%f, ", e, gaus2d_fit[e]);

  // Using the Nelder-Mead algorithm
  fprintf (stderr, "\n\n\t...Using Nelder-Mead : \n\t");
  fit2dgaus_neldermead (c, gaus2d_fit);
  make_gaus2d_inspect (c, gaus2d_fit, "fitexamples/0649nm");
  for (int e = 0; e < 6; e++)
    fprintf (stderr, "[%d]=%f, ", e, gaus2d_fit[e]);

  // Load data from file 
  string file2 ("fitexamples/9568_column.ascii");
  gsl_matrix *cold = read_gsl_matrix_ASCII (file2);

  // Do a 2D Fermi fit on the loaded image
  fprintf (stderr, "\n\nFitting cold image with 2D Gaus and 2D Fermi:\n");
  double fermi2d_fit[7] = { 120.7, 4.59, 49.45, 89.64, 111.3, 203.3, -0.137 };
  double gaus2dfit[6] = { 129.65, 31.78, 111.3, 58.64, 203.3, -0.6 };
  //double fermi2d_fit[7] = { 121., 4.0, 88.67, 48.35, 203.0, 111.4, -0.35 };

  // Using the Nelder-Mead algorithm
  fprintf (stderr, "\n\n\t...First using 2D Gaussian Levenberg-Marquardt : \n\t");
  fit2dgaus (cold, gaus2dfit);
  fprintf (stderr, "\n\n\t...Then using 2D Fermi Nelder-Mead : \n\t");
  fit2dfermi_neldermead (cold, fermi2d_fit);

  make_fermi2d_inspect (cold, fermi2d_fit, "fitexamples/9568nm");
  make_fermi2d_gaus2d_inspect (cold, fermi2d_fit, gaus2dfit,
			       "fitexamples/9568nm");
  for (int e = 0; e < 6; e++)
    fprintf (stderr, "[%d]=%f, ", e, fermi2d_fit[e]);


  fprintf (stderr, "\n\n");
  cout << "testing ..." << endl;

  return 0;
}
