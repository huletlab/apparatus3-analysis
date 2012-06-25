
#include "utils/utils.h"

bool VERBOSE;

/* uncommment the main to test the utilities */
int
main (int argc, char **argv)
{
  VERBOSE = true;

  //Reading a fluorescence image from file
  string file ("../examples/testutils/img.fluor");
  double img[ROW][COL];
  cout << ReadFluorImg (file, img) << endl;
  cout << img[0][0] << endl;

  //Reading a .fits image from file and storing it in valarray
  string file2 ("../examples/testutils/1070atoms.fits");
  valarray < unsigned long >imgdata;
  ReadFitsImg (file2, imgdata);

  //Reading a .fits image from file and storing it in a gsl_matrix
  gsl_matrix *m = ReadFitsImg_gsl_matrix (file2);
  cout << "gsl_matrix   rows: " << m->size1 << ", cols: " << m->size2 << endl;

  //Testing the crop image routine
  string reportcrop ("../examples/testutils/reportCROP.INI");
  gsl_matrix *c = cropImage (reportcrop, m);
  cout << "cropped gsl_matrix   rows: " << c->size1 << ", cols: " << c->
    size2 << endl;
  string outputfile ("../examples/testutils/cropped.ASCII");
  save_gsl_matrix_ASCII (c, outputfile);
  cout << "First few elements of cropped image: " << endl;
  for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 5; j++)
	{
	  cout << gsl_matrix_get (c, i, j) << "\t";
	}
      cout << "... " << endl;
    }
  cout << "..." << endl << endl;
  gsl_matrix_free (c);



  //Getting the total number of counts in an image
  cout << "total counts: " << img_counts (m) << endl;
  cout << "First few elements: " << endl;
  for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 5; j++)
	{
	  cout << gsl_matrix_get (m, i, j) << "\t";
	}
      cout << "... " << endl;
    }
  cout << "..." << endl << endl;
  gsl_matrix_free (m);
  return 0;
}
