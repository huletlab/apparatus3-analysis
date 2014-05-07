/*
 * Project:  Testing fourier transform routines
 *
 * Author:   Pedro M Duarte 2011-02
 * 
 */

#include "fourier.h"
#include "utils/utils.h"


using namespace std;

bool VERBOSE;

int
main (int argc, char **argv)
{

  VERBOSE = false;

  string file ("4504atoms.fits");
  gsl_matrix *img = ReadFitsImg_gsl_matrix (file);

  cout << "testing ..." << endl;
  gsl_matrix *FT = fft2d (img);

  cout << "testing ..." << endl;

  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string ft_path ("");
  ft_path += path;
  ft_path += "/";
  ft_path += "fft2d_";
  ft_path += "4504";
  ft_path += ".TIFF";
  toTiffImage (FT, ft_path);

  cout << "testing ..." << endl;
  gsl_matrix_free (FT);

  cout << "testing ..." << endl;

  return 0;
}
