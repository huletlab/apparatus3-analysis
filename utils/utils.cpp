/*
 * Project:  Utilities for data anlysis. Mostly functions to read data or 
 *           information from various data formats.
 *
 * File:     utils.cpp
 *
 * Author:   Pedro M Duarte 2009-04
 * 
 */

#include "utils/utils.h"

/********** FILE UTILITIES **********/

extern bool VERBOSE;

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

//Reads a fluorescence image from file
bool
ReadFluorImg (string & datafile, double img[ROW][COL])
{
  bool errflg = false;
  ifstream in (datafile.c_str ());
  cout << datafile << " is good: " << in.good () << endl;
  for (int i = 0; i < ROW; i++)
    {
      for (int j = 0; j < COL; j++)
	{
	  if (!in.good ())
	    errflg = true;
	  else
	    in >> img[i][j];
	}
    }
  in.close ();
  return errflg;
}

gsl_matrix *
ReadFluorImg_gsl_matrix (string & datafile)
{
  bool errflg = false;
  ifstream in (datafile.c_str ());
  //cout << datafile << " is good: " << in.good () << endl;

  gsl_matrix *imgdata = gsl_matrix_alloc (1034, 779);	//rows,columns
  double elem;

  for (int i = 0; i < 1034; i++)
    {
      for (int j = 0; j < 779; j++)
	{
	  if (!in.good ())
	    errflg = true;
	  else
	    in >> elem;
	  gsl_matrix_set (imgdata, i, 779 - 1 - j, elem);
//        cout << "i=" << i << ", j=" << j << ",  m(i,j)=" << elem << endl;
//        cin.get ();
	}
    }
  in.close ();
  return imgdata;
}

//Read a .fits image from file and store it in a valarray
bool
ReadFitsImg (string & datafile, valarray < unsigned long >&imgdata)
{
  CCfits::FITS f (datafile);
  CCfits::PHDU & image = f.pHDU ();
  image.readAllKeys ();
  image.read (imgdata);
  //Display header information
  if (VERBOSE)
    cout << image << endl;
  // If you want to check the dimensions of the image:
  //  if ( ix != dimx || iy != dimy ) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

//Read a .fits image from file and store it in a gsl_matrix
gsl_matrix *
ReadFitsImg_gsl_matrix (string & datafile)
{
  CCfits::FITS f (datafile);
  CCfits::PHDU & image = f.pHDU ();
  image.readAllKeys ();
  valarray < unsigned long >contents;
  image.read (contents);
  //Display header information
  if (VERBOSE && false)
    cout << image << endl;
  int ix = (int) image.axis (0);
  int iy = (int) image.axis (1);

  // If you want to check the dimensions of the image:
  //  if ( ix != dimx || iy != dimy ) return EXIT_FAILURE;


  gsl_matrix *imgdata = gsl_matrix_alloc (ix, iy);	//rows,columns


  for (int i = 0; i < ix; i++)
    {
      for (int j = 0; j < iy; j++)
	{
	  gsl_matrix_set (imgdata, i, j, contents[j * ix + i]);
	}
    }
  return imgdata;
}

bool
save_gsl_matrix_ASCII (gsl_matrix * m, string & file)
{
  FILE *fileout;
  fileout = fopen (file.c_str (), "w+");
  if (VERBOSE)
    {
      cout << "\nsaving matrix to ASCII: " << file << endl;
    }
  for (unsigned int j = 0; j < m->size2; j++)
    {
      for (unsigned int i = 0; i < m->size1; i++)
	{
	  fprintf (fileout, "%.5f\t", gsl_matrix_get (m, i, j));
	}
      fprintf (fileout, "\n");
    }
  fclose (fileout);
  return EXIT_SUCCESS;
}


void
toTiffImage (gsl_matrix * m, string & filename, bool invert)
{
  TIFF *tif = TIFFOpen (filename.c_str (), "w");

  unsigned int xWidth = m->size1, yWidth = m->size2;
  TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, xWidth);
  TIFFSetField (tif, TIFFTAG_IMAGELENGTH, yWidth);
  TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 8);
  //TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_LEFTTOP);
  //TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  char *buf = new char[xWidth];

  double min = 1e10, max = 0.0, elem;
  for (unsigned int j = 0; j < yWidth; j++)
    {
      for (unsigned int i = 0; i < xWidth; i++)
	{
	  elem = gsl_matrix_get (m, i, j);
	  if (elem < min)
	    min = elem;
	  if (elem > max)
	    max = elem;
	}
    }

  for (unsigned int j = 0; j < yWidth; j++)
    {
      for (unsigned int i = 0; i < xWidth; i++)
	{
	  elem = gsl_matrix_get (m, i, yWidth - 1 - j);

	  // Here the value for the pixel is normalized to be between 0 and 256
	  // choose one or the other to get white to black or black to white gradients
	  if (invert == true)
	    {
	      buf[i] = char (255. * (1. - (elem - min) / (max - min)));
	    }
	  else
	    {
	      buf[i] = char (255. * ((elem - min) / (max - min)));
	    }
	}
      TIFFWriteScanline (tif, buf, j);
    }

  TIFFClose (tif);

}

//This function is used to created a column data file, where 
//each column is one of the vectors given
void
to_dat_file (gsl_vector * vecs[], unsigned int N, string shot, string datfile)
{

  unsigned int minsize = 16000;
  for (int i = 0; i < (int) N; i++)
    {
      if (vecs[i]->size < minsize)
	minsize = vecs[i]->size;
    }

  unsigned int maxsize = 0;
  for (int i = 0; i < (int) N; i++)
    {
      if (vecs[i]->size > maxsize)
	maxsize = vecs[i]->size;
    }

  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string fname (path);
  fname += "/";
  fname += shot;
  fname += "_";
  fname += datfile;

  FILE *dat;
  dat = fopen (fname.c_str (), "w+");

  for (unsigned int index = 0; index < maxsize; index++)
    {
      for (int i = 0; i < (int) N; i++)
	{
	  if (index >= vecs[i]->size)
	    fprintf (dat, "nan\t");
	  //fprintf (dat, "%e\t", 0.0);
	  else
	    fprintf (dat, "%e\t", gsl_vector_get (vecs[i], index));
	}
      fprintf (dat, "\n");
    }
  fclose (dat);
  return;
}




/********** IMAGE MANIPULATION UTILITIES **********/

void
getmaxRowCol (gsl_matrix * m, gsl_vector * max_row, gsl_vector * max_col)
{

  unsigned int s1 = m->size1;
  unsigned int s2 = m->size2;

  double max = 0.;
  unsigned int i_max = 0;
  unsigned int j_max = 0;
  findpeak (m, &i_max, &j_max, &max);

  max_row = gsl_vector_alloc (s2);
  for (unsigned int j = 0; j < s2; j++)
    {
      gsl_vector_set (max_row, j, gsl_matrix_get (m, i_max, j));
    }

  max_col = gsl_vector_alloc (s1);
  for (unsigned int i = 0; i < s1; i++)
    {
      gsl_vector_set (max_col, i, gsl_matrix_get (m, i, j_max));
    }

  return;
}

void
findpeak (gsl_matrix * m, unsigned int *i_max_ptr, unsigned int *j_max_ptr,
	  double *max_ptr)
{
  unsigned int ravg = 0;
  findpeak_running_avg (m, i_max_ptr, j_max_ptr, max_ptr, ravg);
  return;
}

void
findpeak_running_avg (gsl_matrix * m, unsigned int *i_max_ptr,
		      unsigned int *j_max_ptr, double *max_ptr,
		      unsigned int ravg)
{

  unsigned int s1 = m->size1;
  unsigned int s2 = m->size2;

  double max = -1.e4;
  unsigned int i_max = 0;
  unsigned int j_max = 0;
  double elem;

  for (unsigned int i = 0 + ravg; i < s1 - ravg; i++)
    {
      for (unsigned int j = 0 + ravg; j < s2 - ravg; j++)
	{
	  if (ravg == 0)
	    {
	      elem = gsl_matrix_get (m, i, j);
	    }
	  else
	    {
	      elem = 0.;
	      for (unsigned int ii = i - ravg; ii <= i + ravg; i++)
		{
		  for (unsigned int jj = j - ravg; jj <= j + ravg; j++)
		    {
		      cout << ii << "\t";

		      elem +=
			gsl_matrix_get (m, ii, jj) / pow (2. * ravg + 1., 2.);
		}}
	    }
	  if (elem > max)
	    {
	      max = elem;
	      i_max = i;
	      j_max = j;

	    }
	}
    }

  if (max == 0. || i_max == 0 || j_max == 0)
    {
      cout << "error finding peak: could not find peak" << endl;
      return;
    }
  *i_max_ptr = i_max;
  *j_max_ptr = j_max;
  *max_ptr = max;
  if (VERBOSE)
    {
      cout << endl << "\tPeak found at ( " << i_max << ", " << j_max <<
	" ) = " << max;
    }
  return;
}

unsigned int
to_uint (double x)
{
  double out, temp;
  if (modf (x, &temp) >= 0.5)
    out = x >= 0 ? ceil (x) : floor (x);
  else
    out = x < 0 ? ceil (x) : floor (x);
  // cout << "rounded to " << out << endl;
  return (unsigned int) out;
}







void
findmoments (gsl_matrix * m, unsigned int *ci, unsigned int *cj, double *peak,
	     unsigned int *wi1e, unsigned int *wj1e)
{

  findcenter (m, ci, cj, peak);

  //sometimes a negative atom number shows up in the column density
  //this affects the center of mass calculation, so whenever a pixel
  //with a negative number of atoms is found it is taken as zero 
  //for the center of mass estimate
  double i0 = (double) *ci;
  double j0 = (double) *cj;

  double masstotal = 0., mass = 0.;
  double wi1e_ = 0., wj1e_ = 0.;

  for (unsigned int i = 0; i < m->size1; i++)
    for (unsigned int j = 0; j < m->size2; j++)
      {
	mass = gsl_matrix_get (m, i, j);
	mass = mass > 0 ? mass : 0.;
	masstotal += mass;
	wi1e_ += mass * (i - i0) * (i - i0);
	wj1e_ += mass * (j - j0) * (j - j0);
      }
  wi1e_ = sqrt (2. * wi1e_ / masstotal);
  wj1e_ = sqrt (2. * wj1e_ / masstotal);

/*
//// SUM ALONG I AND J FIRST
  gsl_vector *isum = gsl_vector_alloc (m->size1);
  gsl_vector *jsum = gsl_vector_alloc (m->size1);

  double sum;

  for (unsigned int i = 0; i < m->size1; i++)
    {
      sum = 0.;
      for (unsigned int j = 0; j < m->size2; j++)
	{
	  sum = sum + gsl_matrix_get (m, i, j);
	}
      gsl_vector_set (isum, i, sum);
    }

  for (unsigned int j = 0; j < m->size2; j++)
    {
      sum = 0.;
      for (unsigned int i = 0; i < m->size1; i++)
	{
	  sum = sum + gsl_matrix_get (m, i, j);
	}
      gsl_vector_set (jsum, j, sum);
    }

  double masstotal = 0., mass = 0.;
  double wi1e_ = 0., wj1e_ = 0.;

  for (unsigned int i = 0; i < m->size1; i++)
    {
      mass = gsl_vector_get (isum, i);
      mass = mass > 0 ? mass : 0.;
      masstotal += mass;
      wi1e_ += mass * (i - i0) * (i - i0);
    }
  wi1e_ = sqrt (2. * wi1e_ / masstotal);


  masstotal = 0.;
  mass = 0.;
  for (unsigned int j = 0; j < m->size2; j++)
    {
      mass = gsl_vector_get (jsum, j);
      mass = mass > 0 ? mass : 0.;
      masstotal += mass;
      wj1e_ += mass * (j - j0) * (j - j0);
    }
  wj1e_ = sqrt (2. * wj1e_ / masstotal);

////
*/

  if (VERBOSE && false)
    {
      printf ("Center is at (%.1f,%.1f)\n", i0, j0);
      cout << "Moments results:  ci = " << i0 << ",  cj = "
	<< j0 << ",  wi1e = " << wi1e_ << ",  wj1e = " << wj1e_ << endl <<
	endl;;
    }

  *wi1e = (unsigned int) floor (wi1e_);
  *wj1e = (unsigned int) floor (wj1e_);

  return;
}


void
findcenter (gsl_matrix * m, unsigned int *i_max_ptr, unsigned int *j_max_ptr,
	    double *max_ptr)
{
  double masstotal = 0, mass = 0;

  double ci = 0.;
  double cj = 0.;

  //sometimes a negative atom number shows up in the column density
  //this affects the center of mass calculation, so whenever a pixel
  //with a negative number of atoms is found it is taken as zero 
  //for the center of mass estimate

  for (unsigned int i = 0; i < m->size1; i++)
    for (unsigned int j = 0; j < m->size2; j++)
      {
	mass = gsl_matrix_get (m, i, j);
	mass = mass > 0 ? mass : 0.;
	masstotal += mass;
	ci += mass * i;
	cj += mass * j;
      }
  ci = ci / masstotal;
  cj = cj / masstotal;

  if (VERBOSE)
    {
      cout << endl << "Center of mass results:  ci = " << ci << ",  cj = "
	<< cj << endl;
    }

  unsigned int i_max = 0;
  unsigned int j_max = 0;
  double max = -1.e4;


  i_max = to_uint (ci);
  j_max = to_uint (cj);

  if (VERBOSE)
    {
      cout << "Conversion to unsigned int:  ci = " << i_max << ",  cj = "
	<< j_max << endl;
    }

  max = gsl_matrix_get (m, i_max, j_max);

  if (i_max == 0 || j_max == 0)
    {
      cout << "error finding center: could not find center" << endl;
      return;
    }
  if (max <= 0.)
    {
      if (VERBOSE)
	cout << endl <<
	  "Warning: center pixel has negative atom number, max overridden" <<
	  endl;
      max = 10.;
    }
  *i_max_ptr = i_max;
  *j_max_ptr = j_max;
  *max_ptr = max;
  if (VERBOSE)
    {
      cout << "\tCenter found at ( " << i_max << ", " << j_max <<
	" ) = " << max << endl;;
    }
  return;
}


void
findFWHM (gsl_matrix * m, unsigned int *FWHM_i, unsigned int *FWHM_j)
{
  unsigned int s1 = m->size1;
  unsigned int s2 = m->size2;

  double max = 0.;
  unsigned int i_max = 0;
  unsigned int j_max = 0;
  //findpeak (m, &i_max, &j_max, &max);
  findcenter (m, &i_max, &j_max, &max);

  double elem;

  unsigned int i_FWHM_plus = i_max;
  while (i_FWHM_plus < s1)
    {
      elem = gsl_matrix_get (m, i_FWHM_plus, j_max);
      if (elem < max / 2.)
	break;
      i_FWHM_plus++;
    }
  unsigned int i_FWHM_minus = i_max;
  while (i_FWHM_minus > 0)
    {
      elem = gsl_matrix_get (m, i_FWHM_minus, j_max);
      if (elem < max / 2.)
	break;
      i_FWHM_minus--;
    }

  unsigned int j_FWHM_plus = j_max;
  while (j_FWHM_plus < s2)
    {
      elem = gsl_matrix_get (m, i_max, j_FWHM_plus);
      if (elem < max / 2.)
	break;
      j_FWHM_plus++;
    }
  unsigned int j_FWHM_minus = j_max;
  while (j_FWHM_minus > 0)
    {
      elem = gsl_matrix_get (m, i_max, j_FWHM_minus);
      if (elem < max / 2.)
	break;
      j_FWHM_minus--;
    }

  *FWHM_i = i_FWHM_plus - i_FWHM_minus;
  *FWHM_j = j_FWHM_plus - j_FWHM_minus;

  if (*FWHM_i == 0)
    {
      if (VERBOSE)
	cout << "Warning: FWHM_i overridden to be > 0" << endl;
      *FWHM_i = 10;
    }
  if (*FWHM_j == 0)
    {
      if (VERBOSE)
	cout << "Warning: FWHM_i overridden to be > 0" << endl;
      *FWHM_j = 10;
    }

  if (VERBOSE)
    {
      cout << endl << "\tFWHM_i = " << *FWHM_i << ",  FWHM_j = " << *FWHM_j <<
	endl;
    }
  return;
}

gsl_matrix *
mask (gsl_matrix * m, double factor)
{
  double min = 1e10, max = 0.0, elem;
  for (unsigned int i = 0; i < m->size1; i++)
    {
      for (unsigned int j = 0; j < m->size2; j++)
	{
	  elem = gsl_matrix_get (m, i, j);
	  if (elem < min)
	    min = elem;
	  if (elem > max)
	    max = elem;
	}
    }

  gsl_matrix *masked = gsl_matrix_alloc (m->size1, m->size2);
  for (unsigned int i = 0; i < m->size1; i++)
    {
      for (unsigned int j = 0; j < m->size2; j++)
	{
	  if (gsl_matrix_get (m, i, j) > min + factor * (max - min))
	    gsl_matrix_set (masked, i, j, 1.0);
	  else
	    gsl_matrix_set (masked, i, j, 0.);
	}
    }
  return masked;
}

gsl_matrix *
smooth (gsl_matrix * raw, unsigned int bins)
{

  if (bins < 2)
    return raw;
  unsigned int s1 = raw->size1;
  unsigned int s2 = raw->size2;
  gsl_matrix *smoothed = gsl_matrix_alloc (s1, s2);

  double avg;
  int count;
  unsigned int lmin, lmax, mmin, mmax;
  for (unsigned int i = 0; i < raw->size1; i++)
    {
      for (unsigned int j = 0; j < raw->size2; j++)
	{
	  //printf ("(%d,%d) -->", i, j);
	  count = 0;
	  avg = 0.;

	  lmin = (unsigned int) max ((double) i - (double) bins, 0.);
	  lmax =
	    (unsigned int) min ((double) i + (double) bins, (double) s1 - 1.);

	  mmin = (unsigned int) max ((double) j - (double) bins, 0.);
	  mmax =
	    (unsigned int) min ((double) j + (double) bins, (double) s2 - 1.);

	  //printf ("(%d,%d) to (%d,%d)\n", lmin, lmax, mmin, mmax);
	  for (unsigned int l = lmin; l <= lmax; l++)
	    {
	      for (unsigned int m = mmin; m <= mmax; m++)
		{
		  avg += gsl_matrix_get (raw, l, m);
		  count++;
		}
	    }
	  avg = avg / count;
	  gsl_matrix_set (smoothed, i, j, avg);
	}
    }
  return smoothed;
}

gsl_matrix *
autocropImage (gsl_matrix * raw, double nFWHM)
{
  unsigned int s1 = raw->size1;
  unsigned int s2 = raw->size2;

  double max = 0.;
  unsigned int i_max = 0;
  unsigned int j_max = 0;
  findpeak (raw, &i_max, &j_max, &max);

  if (max == 0. || i_max == 0 || j_max == 0)
    {
      cout << "error finding peak: could not find peak" << endl;
      return NULL;
    }

  unsigned int FWHM_i, FWHM_j;
  findFWHM (raw, &FWHM_i, &FWHM_j);

  unsigned int i0 = (unsigned int) floor (i_max - FWHM_i * nFWHM / 2.);
  unsigned int i1 = (unsigned int) floor (i_max + FWHM_i * nFWHM / 2.);

  unsigned int j0 = (unsigned int) floor (j_max - FWHM_j * nFWHM / 2.);
  unsigned int j1 = (unsigned int) floor (j_max + FWHM_j * nFWHM / 2.);

  if (i0 < 0 || i0 > s1 || i1 < 0 || i1 > s1 || i0 > i1)
    {
      cout << "error finding peak: i is out of bounds" << endl;
      return NULL;
    }
  if (j0 < 0 || j0 > s2 || j1 < 0 || j1 > s2 || j0 > j1)
    {
      cout << "error finding peak: j is out of bounds" << endl;
      return NULL;
    }

  gsl_matrix *cropped = gsl_matrix_alloc (i1 - i0, j1 - j0);
  for (unsigned int i = i0; i < i1; i++)
    {
      for (unsigned int j = j0; j < j1; j++)
	{
	  gsl_matrix_set (cropped, i - i0, j - j0,
			  gsl_matrix_get (raw, i, j));
    }}

  return cropped;
}

gsl_matrix *
cropImage_ROI (unsigned int roi[], gsl_matrix * raw)
{


  int ax0pos = roi[0];
  int ax1pos = roi[1];
  int ax0size = roi[2];
  int ax1size = roi[3];

  if (VERBOSE)
    {
      cout << endl;
      cout << "Cropping image (" << raw->size1 << ", " << raw->size2 << ")";
      cout << " ==> ( " << ax0pos << ":" << ax0pos << "+" << ax0size << ", "
	<< ax1pos << ":" << ax1pos << "+" << ax1size << " )";
    }

  gsl_matrix *cropped = gsl_matrix_alloc (ax0size, ax1size);


  for (int i = ax0pos; i < ax0pos + ax0size; i++)
    {
      for (int j = ax1pos; j < ax1pos + ax1size; j++)
	{
	  gsl_matrix_set (cropped, i - ax0pos, j - ax1pos,
			  gsl_matrix_get (raw, i, j));
    }}

  return cropped;
}

gsl_matrix *
cropImage (string & reportfile, gsl_matrix * raw)
{


  int ax0pos = (int) getINI_num (reportfile, "ROI", "ax0pos");
  int ax1pos = (int) getINI_num (reportfile, "ROI", "ax1pos");
  int ax0size = (int) getINI_num (reportfile, "ROI", "ax0size");
  int ax1size = (int) getINI_num (reportfile, "ROI", "ax1size");

  if (VERBOSE)
    {
      cout << endl;
      cout << "Cropping image (" << raw->size1 << ", " << raw->size2 << ")";
      cout << " ==> ( " << ax0pos << ":" << ax0pos << "+" << ax0size << ", "
	<< ax1pos << ":" << ax1pos << "+" << ax1size << " )";
    }

  gsl_matrix *cropped = gsl_matrix_alloc (ax0size, ax1size);


  for (int i = ax0pos; i < ax0pos + ax0size; i++)
    {
      for (int j = ax1pos; j < ax1pos + ax1size; j++)
	{
	  gsl_matrix_set (cropped, i - ax0pos, j - ax1pos,
			  gsl_matrix_get (raw, i, j));
    }}

  return cropped;
}

gsl_matrix *
subtract (gsl_matrix * m1, gsl_matrix * m2)
{
  if ((m1->size1 != m2->size1) || (m1->size2 != m2->size2))
    {
      cout <<
	"Attempted subraction of two images with different dimensions." <<
	endl;
      cout << "Program will exit " << endl;
    }
  //cout << m1->size1 << " , " << m1->size2 << endl;
  gsl_matrix *s = gsl_matrix_alloc (m1->size1, m1->size2);
  for (unsigned int i = 0; i < m1->size1; i++)
    {
      for (unsigned int j = 0; j < m1->size2; j++)
	{
	  double el = gsl_matrix_get (m1, i, j) - gsl_matrix_get (m2, i, j);
	  gsl_matrix_set (s, i, j, el);
	}
    }
  return s;
}

unsigned int
coerce_matrix_index (unsigned int i, unsigned int size)
{
  if (VERBOSE)
    {
      cout << "--- coercing gsl_matrix index: " << i << " into " << size <<
	endl;
    }
  if (i >= size)
    return size - 1;
  if (i < 0)
    return 0;
  return i;
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
