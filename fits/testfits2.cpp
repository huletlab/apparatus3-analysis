/*
 * Project:  Testing fitting routines
 *
 * Author:   Pedro M Duarte 2011-02
 * 
 */

#include "funcs/funcs.h"
#include "utils/utils.h"
#include "fits/fits.h"

#include <math.h>

using namespace std;

bool VERBOSE;
bool DEBUG_FUNCS;
bool DEBUG_FITS;

int
main (int argc, char **argv)
{

  VERBOSE = true;

  // Load data from file 
  string file2 ("fitexamples/0102_column.ascii");
  gsl_matrix *cold = read_gsl_matrix_ASCII (file2);

  // Do a 2D Fermi fit on the loaded image
  fprintf (stderr, "\n\nFitting cold image with 2D Gaus and 2D Fermi:\n");

  double fermi2d_fit[7] = { 69.55, -4.00, 18.72, 48.23, 65.76, 168.64, 0.0038};
  double gaus2dfit[6] = { 65.75, 18.77, 168.61, 48.37, 69.11, -0.003 };

  // Using the Nelder-Mead algorithm
  fprintf (stderr,
	   "\n\n\t...First using 2D Gaussian Levenberg-Marquardt : \n\t");
  fit2dgaus (cold, gaus2dfit);
  fprintf (stderr, "\n\n\t...Then using 2D Fermi Nelder-Mead : \n\t");
  fit2dfermi_neldermead (cold, fermi2d_fit);

  make_fermi2d_inspect (cold, fermi2d_fit, "9568nm", "fitexamples/9568nm");
  make_fermi2d_gaus2d_inspect (cold, fermi2d_fit, gaus2dfit, "9568nm",
			       "fitexamples/9568nm");

  //Try out the High Level 2D Gaus fit
  double gaus2dfit_high_level[6];
  double gaus2dfit_high_level_err[6];
  string prefix ("fitexamples/high_level");
  Fit2DGaus_High_Level (cold, gaus2dfit_high_level, gaus2dfit_high_level_err,
			prefix);

  fprintf (stderr, "\n\nResults of first 2D Gaus Fit\n");
  for (int e = 0; e < 6; e++)
    fprintf (stderr, "[%d]=%f, ", e, gaus2dfit[e]);

  fprintf (stderr, "\n\nResults of second  2D Gaus Fit\n");
  for (int e = 0; e < 6; e++)
    fprintf (stderr, "[%d]=%f, ", e, gaus2dfit_high_level[e]);

  fprintf (stderr, "\n\nResults of 2D Fermi Fit\n");
  for (int e = 0; e < 7; e++)
    fprintf (stderr, "[%d]=%f, ", e, fermi2d_fit[e]);


  fprintf (stderr, "\n\n");
  cout << "testing ..." << endl;

  return 0;
}
