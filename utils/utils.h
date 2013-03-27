#ifndef UTILS_HEADER_FILE
#define UTILS_HEADER_FILE

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


/* 
 * Functions in fileutils.cpp that will be exported:
 *
 */
string makepath (const char *base, string prefix, const char *identifier);
int makeShotPaths_Basler (char *shot, string & shotnum, string & report,
			  string & atoms);
int makeShotPaths (char *shot, string & shotnum, string & report,
		   string & atoms, string & noatoms, string & atomsref,
		   string & noatomsref);

int NLines (string datafile);



/* 
 * Functions in readwrite.cpp that will be exported:
 *
 */
#define ROW 779
#define COL 1034

bool ReadFluorImg (string & datafile, double img[ROW][COL]);
gsl_matrix *ReadFluorImg_gsl_matrix (string & datafile);
bool ReadFitsImg (string & datafile, valarray < unsigned long >&imgdata);
gsl_matrix *ReadZeros_gsl_matrix (const gsl_matrix * shape);
gsl_matrix *ReadFitsImg_gsl_matrix (string & datafile);

bool save_gsl_matrix_ASCII (gsl_matrix * m, string & file);
gsl_matrix *read_gsl_matrix_ASCII (string & file);
void toTiffImage (gsl_matrix * m, string & filename, bool invert = false);

void to_dat_file (gsl_vector * vecs[], unsigned int N, string prefix,
		  string sufix);

void
to_dat_file_2 (gsl_vector * one, gsl_vector * two, string prefix,
	       string sufix);


/* 
 * Functions in matrices.cpp that will be exported:
 *
 */
void getmaxRowCol (gsl_matrix * m, gsl_vector * max_row,
		   gsl_vector * max_col);
void findpeak (gsl_matrix * m, unsigned int *i_max_ptr,
	       unsigned int *j_max_ptr, double *max_ptr, bool is_pos = true);
void findpeak_running_avg (gsl_matrix * m, unsigned int *i_max_ptr,
			   unsigned int *j_max_ptr, double *max_ptr,
			   unsigned int ravg, bool is_pos = true);
void findmoments (gsl_matrix * m, unsigned int *ci, unsigned int *cj,
		  double *peak, unsigned int *wi1e, unsigned int *wj1e);
void findcenter (gsl_matrix * m, unsigned int *i_max_ptr,
		 unsigned int *j_max_ptr, double *max_ptr);
void findFWHM (gsl_matrix * m, unsigned int *FWHM_i, unsigned int *FWHM_j);

gsl_matrix *mask (gsl_matrix * m, double factor = 5);
gsl_matrix *smooth (gsl_matrix * raw, unsigned int bins);
gsl_matrix *subtract (gsl_matrix * m1, gsl_matrix * m2);

/*void Gaus2DGuess (gsl_matrix * m, double *guess, string prefix,
		  bool save_matrices);*/

unsigned int coerce_matrix_index (unsigned int i, unsigned int size);

double img_counts (gsl_matrix * m);
double img_peak (gsl_matrix * m, double *pos);



/* 
 * Functions in crop.cpp that will be exported:
 *
 */
gsl_matrix *autocropImage (gsl_matrix * raw, double nFWHM);
gsl_matrix *cropImage_ROI (unsigned int roi[], gsl_matrix * raw);
gsl_matrix *cropImage (string & reportfile, gsl_matrix * raw);


#endif
