/*
 * Project:  This file contain functions to manipulate matrices
 *            
 *
 * Author:   Pedro M Duarte 2012-06
 * 
 */


extern bool VERBOSE;

/***** Function prototypes *****/


void getmaxRowCol(gsl_matrix *m, gsl_vector * max_row, gsl_vector * max_col); 
void findpeak( gsl_matrix *m, unsigned int * i_max_ptr, unsigned int * j_max_ptr, double * max_ptr); 
void findpeak_running_avg( gsl_matrix *m, unsigned int * i_max_ptr, unsigned int * j_max_ptr, double * max_ptr, unsigned int ravg);
void findmoments(gsl_matrix *m, unsigned int *ci, unsigned int *cj,  double *peak, unsigned int *wi1e, unsigned int *wj1e);
void findcenter( gsl_matrix *m, unsigned int * i_max_ptr, unsigned int * j_max_ptr, double * max_ptr); 
void findFWHM ( gsl_matrix * m, unsigned int * FWHM_i, unsigned int * FWHM_j);

gsl_matrix *mask ( gsl_matrix *m, double factor=5);
gsl_matrix * smooth(gsl_matrix * raw, unsigned int bins);  
gsl_matrix * subtract( gsl_matrix* m1, gsl_matrix* m2);

unsigned int coerce_matrix_index( unsigned int i, unsigned int size);

 
/***** Function code *****/



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

