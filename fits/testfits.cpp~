/*
 * Project:  Testing fitting routines
 *
 * Author:   Pedro M Duarte 2011-02
 * 
 */

#include "funcs/funcs.h"
#include "utils/utils.h"
#include "fits/fits.h"
#include <gsl/gsl_sf_erf.h>
#include <math.h>

using namespace std;

bool VERBOSE;
bool DEBUG_FUNCS;
bool DEBUG_FITS;

int
main (int argc, char **argv)
{

  VERBOSE = true;
/*
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
    fprintf (stderr, "[%d]=%f, ", e, gaus2d_fit[e]); */


  // Load data from file 
  string file2 ("fitexamples/9568_column.ascii");
  gsl_matrix *cold = read_gsl_matrix_ASCII (file2);

  // Do a 2D Fermi fit on the loaded image
  fprintf (stderr, "\n\nFitting cold image with 2D Gaus and 2D Fermi:\n");
  double fermi2d_fit[7] = { 120.7, 4.59, 49.45, 89.64, 111.3, 203.3, -0.137 };
  double gaus2dfit[6] = { 111.3, 31.78, 203.3, 58.64, 129.8, -0.6 };
  double mottgaus2dfit[6] = { 111.3, 31.78, 203.3, 58.64, 129.8, 0 };
  //double fermi2d_fit[7] = { 121., 4.0, 88.67, 48.35, 203.0, 111.4, -0.35 };

  // Using the Nelder-Mead algorithm
  fprintf (stderr,
	   "\n\n\t...First using 2D Gaussian Levenberg-Marquardt : \n\t");
  fit2dgaus (cold, gaus2dfit);

  fprintf (stderr, "\n\n\t...Using 2D MOTT Gaussian simplex : \n\t");
  fit2dmottgaus_neldermead (cold, mottgaus2dfit);
  printf
    ("Gaus2d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e ",
     gaus2dfit[0],
     gaus2dfit[1], gaus2dfit[2], gaus2dfit[3], gaus2dfit[4], gaus2dfit[5]);
  /*
     fprintf (stderr, "\n\n\t...Then using 2D Fermi Nelder-Mead : \n\t");

     fit2dfermi_neldermead (cold, fermi2d_fit);

     make_fermi2d_inspect (cold, fermi2d_fit, "fitexamples/9568nm");
     make_fermi2d_gaus2d_inspect (cold, fermi2d_fit, gaus2dfit,
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
   */

  fprintf (stderr, "\n\n");
  cout << "testing ..." << endl;
  ofstream outputfile;
  outputfile.open ("testdata.txt");
  gsl_vector *g2d = gsl_vector_alloc (6);
  gsl_vector_set (g2d, 0, 0);
  gsl_vector_set (g2d, 1, gaus2dfit[1]);
  gsl_vector_set (g2d, 2, 0);
  gsl_vector_set (g2d, 3, gaus2dfit[3]);
  gsl_vector_set (g2d, 4, gaus2dfit[4]);
  gsl_vector_set (g2d, 5, gaus2dfit[5]);
  gsl_vector *mg2d = gsl_vector_alloc (6);
  gsl_vector_set (mg2d, 0, 0);
  gsl_vector_set (mg2d, 1, mottgaus2dfit[1]);
  gsl_vector_set (mg2d, 2, 0);
  gsl_vector_set (mg2d, 3, mottgaus2dfit[3]);
  gsl_vector_set (mg2d, 4, mottgaus2dfit[4]);
  gsl_vector_set (mg2d, 5, abs (mottgaus2dfit[5]));
  cout << gsl_vector_get (mg2d, 5) << endl;
  for (int i = 0; i < 3 * mottgaus2dfit[1]; i++)
    {
      double x = i * 1.0;
      outputfile << x << "\t" << gaus2d_model (x, 0, g2d) << "\t" <<
	mottgaus2d_model (x, 0, mg2d) << "\t" <<
	gsl_matrix_get ((gsl_matrix *) cold, int (x + gaus2dfit[0]),
			int (gaus2dfit[2])) << endl;
    }
  outputfile.close ();
  double nfit = gaus2dfit[4] * M_PI * gaus2dfit[1] * gaus2dfit[3];
  double r0 = abs (mottgaus2dfit[5]);
  double nfit_mott =
    mottgaus2dfit[4] * M_PI * mottgaus2dfit[1] * mottgaus2dfit[3] * (1 +
								     2 /
								     M_SQRTPI
								     * r0 *
								     exp (-r0
									  *
									  r0)
								     -
								     gsl_sf_erf
								     (r0) +
								     4.0 /
								     3.0 /
								     M_SQRTPI
								     *
								     pow (r0,
									  3) *
								     exp (-r0
									  *
									  r0));
  double peakd = gaus2dfit[4] / (pow (M_PI, 0.5) * gaus2dfit[1]) / pow (1 * 1e-4, 3);	// cm^-3
  double peakd_mott =
    mottgaus2dfit[4] / M_SQRTPI / mottgaus2dfit[1] * exp (-r0 * r0) / pow (1 *
									   1e-4,
									   3);
  cout << "nfit\tnift_mott\tpeakd\tpeakd_mott\n" << nfit << "\t" << nfit_mott
    << "\t" << peakd << "\t" << peakd_mott << endl;

  return 0;
}
