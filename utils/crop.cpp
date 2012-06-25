
/*
 * Project:  This file contain functions to apply crpping procedures to 
 *           matrices 
 *
 * Author:   Pedro M Duarte 2012-06
 * 
 */

#include "utils/utils.h"

extern bool VERBOSE;


/***** Function code *****/


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
