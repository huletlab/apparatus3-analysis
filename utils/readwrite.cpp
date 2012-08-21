/*
 * Project:  This file contain functions to read/write various kinds of
 *           data to various kinds of file formats
 *            
 *
 * Author:   Pedro M Duarte 2012-06
 * 
 */


#include "utils/utils.h"

extern bool VERBOSE;



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
	  if (!in.good () && errflg == false)
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

gsl_matrix *
read_gsl_matrix_ASCII (string & file)
{
  ifstream in (file.c_str ());
  //cout << datafile << " is good: " << in.good () << endl;

  string line;
  int RowNum = 0;
  int ColNum;
  while (getline (in, line, '\n'))
    {
      string entry;
      int EntriesPerRow = 0;
      stringstream ss;
      ss << line;
      while (getline (ss, entry, '\t'))
	EntriesPerRow++;
      if (RowNum == 0)
	ColNum = EntriesPerRow;
      else if (EntriesPerRow != ColNum)
	{
	  cerr << " ERROR Reading matrix ASCII. Program will stop" << endl;
	  exit (1);
	}
      RowNum++;
    }

  if (VERBOSE)
    {
      cout << endl << "Reading gsl_matrix from ASCII file : " << file << endl;
      printf ("\t lines = %d, entries per line = %d\n", RowNum, ColNum);
    }
  in.clear ();
  in.seekg (0, ios::beg);
  gsl_matrix *m = gsl_matrix_alloc (ColNum, RowNum);

  for (unsigned int j = 0; j < m->size2; j++)
    {
      for (unsigned int i = 0; i < m->size1; i++)
	{
	  double elem;
	  in >> elem;
	  gsl_matrix_set (m, i, j, elem);
	}
    }
/*
  unsigned int j = 0;
  while (getline (in, line, '\n'))
    {
      stringstream ss;
      ss << line;
      string entry;
      unsigned int i = 0;
      while (getline (ss, entry, ' '))
	{
	  stringstream entry_ss;
	  double entry_double;
	  char entry_char;
	  entry_ss << entry;
	  entry_ss >> entry_double;

	  if (!entry_ss)
	    {
	      cerr <<
		" ERROR  Reading matrix ASCII (Bad entry).  Program will stop"
		<< endl;
	      exit (1);
	    }
	  else if (entry_ss >> entry_char)
	    {
	      cerr <<
		" ERROR  Reading matrix ASCII (Bad entry).  Program will stop"
		<< endl;
	      exit (1);
	    }
	  else
	    {
	      gsl_matrix_set (m, i, j, entry_double);
	      i++;
	    }
	}
      j++;
    }*/
  return m;
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
to_dat_file (gsl_vector * vecs[], unsigned int N, string prefix, string sufix)
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
  fname += prefix;
  fname += "_";
  fname += sufix;

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

void
to_dat_file_2 (gsl_vector * one, gsl_vector * two, string prefix,
	       string sufix)
{
  gsl_vector *output[2] = { one, two };
  to_dat_file (output, 2, prefix, sufix);
  return;
}
