#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <valarray>

#include "qini/qini.h"
#include "CCfits/config.h"
#include "CCfits/CCfits"
#include "gsl-1.15/gsl/gsl_matrix.h"
#include "gsl-1.15/gsl/gsl_sort_vector.h"
#include "tiff-4.0.0/libtiff/tiffio.h"


#include <math.h>


using namespace std;

int  makeShotPaths_Basler (char *shot, string & shotnum, string & report, string & atoms);
int  makeShotPaths (char *shot, string & shotnum, string & report, string & atoms, string & noatoms, string & atomsref, string & noatomsref);

int  NLines( string datafile);

double img_counts( gsl_matrix *m);
double img_peak ( gsl_matrix *m, double *pos);

double counts2atoms(string &reportfile);


int
makeShotPaths_Basler (char *shot, string & shotnum, string & report,
		      string & atoms)
{

  //fisrt argument is the shot number
  int shot_int = atoi (shot);
  char shot_str[4];
  sprintf (shot_str, "%04d", shot_int);
  shotnum = shot_str;

  //creates reportfile path in current directory
  report = "";
  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string tmp3 ("");
  report += path;
  report += "/report";
  report += shot_str;
  report += ".INI";

  //creates fluor file path in the current directory
  atoms = "";
  atoms += path;
  atoms += "/";
  atoms += shot_str;
  atoms += ".fluor";

  return EXIT_SUCCESS;
}

int
makeShotPaths (char *shot, string & shotnum, string & report,
	       string & atoms, string & noatoms, string & atomsref,
	       string & noatomsref)
{

  //fisrt argument is the shot number
  int shot_int = atoi (shot);
  char shot_str[4];
  sprintf (shot_str, "%04d", shot_int);
  shotnum = shot_str;

  //creates reportfile path in current directory
  report = "";
  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string tmp3 ("");
  report += path;
  report += "/report";
  report += shot_str;
  report += ".INI";

  //creates atoms file path in the current directory
  atoms = "";
  atoms += path;
  atoms += "/";
  atoms += shot_str;
  atoms += "atoms.fits";

  //creates noatoms file path in the current directory
  noatoms = "";
  noatoms += path;
  noatoms += "/";
  noatoms += shot_str;
  noatoms += "noatoms.fits";

  //creates atomsref file path in the current directory
  atomsref = "";
  atomsref += path;
  atomsref += "/";
  atomsref += shot_str;
  atomsref += "atomsref.fits";

  //creates noatomsref file path in the current directory
  noatomsref = "";
  noatomsref += path;
  noatomsref += "/";
  noatomsref += shot_str;
  noatomsref += "noatomsref.fits";

  return EXIT_SUCCESS;
}


//Count number of lines in a text file
int
NLines (string datafile)
{
  ifstream in (datafile.c_str ());
  int size = 0;
  string line;
  while (getline (in, line))
    size++;
  in.close ();
  return size;
}





/********** IMAGE STATISTICS UTILITIES **********/

double
img_counts (gsl_matrix * m)
{
  double c = 0.;
  for (unsigned int i = 0; i < m->size1; i++)
    {
      for (unsigned int j = 0; j < m->size2; j++)
	{
	  c += gsl_matrix_get (m, i, j);
	}
    }
  return c;
}

double
img_peak (gsl_matrix * m, double *pos)
{
  double peak = 0.;
  double pix;
  int imax = 0, jmax = 0;
  for (unsigned int i = 0; i < m->size1; i++)
    {
      for (unsigned int j = 0; j < m->size2; j++)
	{
	  pix = gsl_matrix_get (m, i, j);
//        cout << pix << "\t" << endl;
	  if (pix > peak)
	    {
	      peak = pix;
	      imax = i;
	      jmax = j;
	    }
	}
    }
  pos[0] = imax;
  pos[1] = jmax;
/*  cout << "Peak " << peak << "  at imax=" << imax << ";  jmax=" << jmax <<
  endl; */
  return peak;
}




/********** OTHER UTILITIES **********/

double
counts2atoms (string & reportfile)
{

  double i0 = 6. * getINI_num (reportfile, "ANDOR", "imgpow") / 2.;
  double dettrap = getINI_num (reportfile, "ANDOR", "imgdettrap") / 5.9;
  double detrep = getINI_num (reportfile, "ANDOR", "imgdetrep") / 5.9;
  double powtrap = getINI_num (reportfile, "ANDOR", "imgpowtrap") / 5.9;
  double powrep = getINI_num (reportfile, "ANDOR", "imgpowrep") / 5.9;

  if (powtrap == 0.0 || powrep == 0.0)
    {
      cout <<
	"#We haven't calibrated imaging with a single frequency! Program will exit"
	<< endl;
      exit (EXIT_FAILURE);
    }
  if (dettrap != detrep)
    {
      cout <<
	"#We haven't calibrated imaging with different detuning on trap and repump! Program will exit"
	<< endl;
      exit (EXIT_FAILURE);
    }
  double det = dettrap;

  double texp = getINI_num (reportfile, "ANDOR", "exp");
  double qeff = getINI_num (reportfile, "ANDOR", "quantumeff");
  double sa = getINI_num (reportfile, "ANDOR", "solidangle");
  double g = qeff / sa;

  double c2n =
    (1. / (i0 / (1. + 2. * i0 + 4 * pow (det, 2.)))) * (27.02e-6 / texp) * g;

  return c2n;

}
