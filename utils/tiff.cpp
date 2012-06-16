/* 
	Analysis Program for ODT trap - Time of flight test
	Aurther: Ernie (Acutally is just a modification of Pedro's work)
	Data   : 2011/02/03
*/

#include "utils.h"


struct params
{
  string shotnum;
  string reportfile;
  string atomsfile;
  string noatomsfile;

  bool verbose;
  bool reanalyze;
};

bool VERBOSE;

int
processArgsAnalyze (int argc, char **argv, struct params &p)
{
/*  Read command line arguments */

  VERBOSE = true;

  if (argc == 1)
    {
      cout << endl;
      cout << "  usage:  " << argv[0] << " SHOTNUM [OPTIONS] " << endl;
      cout << endl;
      cout << "  OPTIONS " << endl;
      cout << "\t-v\tverbose" << endl;
      exit (2);
    }

  makeShotPaths (argv[1], p.shotnum, p.reportfile, p.atomsfile,
		 p.noatomsfile);


  p.verbose = false;
  p.reanalyze = false;

  int c;
  while ((c = getopt (argc, argv, "vf")) != -1)
    {
      switch (c)
	{
	case 'v':
	  p.verbose = true;
	case 'f':
	  p.reanalyze = true;
	default:
	  break;
	}
    }

  if (sectionExists (p.reportfile, "CPP") && !p.reanalyze)
    {
      cout << endl;
      cout << " Shot " << p.shotnum << " has already been analyzed." << endl;
      cout << " Use option -f  to force analysis." << endl << endl;
      exit (2);
    }
//  cout << p.shotnum << endl << p.reportfile << endl << p.
//    atomsfile << endl << p.noatomsfile << endl;

  return EXIT_SUCCESS;
}


int
main (int argc, char **argv)
{


  struct params p;
  processArgsAnalyze (argc, argv, p);

  gsl_matrix *atomsNonCrop = ReadFitsImg_gsl_matrix (p.atomsfile);
  gsl_matrix *atoms = cropImage (p.reportfile, atomsNonCrop);
  atomsNonCrop = ReadFitsImg_gsl_matrix (p.noatomsfile);
  gsl_matrix *noatoms = cropImage (p.reportfile, atomsNonCrop);


  stringstream strstream, strstream2, strstream3, strstream4;
  char buffer[50], section[50], para[50];

  if (argc >= 4)
    {
      sprintf (section, argv[3]);
      sprintf (para, argv[4]);
    }
  else
    {
      sprintf (section, "SEQ");
      sprintf (para, "shot");
      cout << section;
    }

  sprintf (buffer, "%05d", int (getINI_num (p.reportfile, section, para)));
  strstream << para << "_" << buffer << "_atom.tiff";
  strstream2 << para << "_" << buffer << "_noatom.tiff";
  string filename, filename2;
  strstream >> filename;
  strstream2 >> filename2;

  toTiffImage (atoms, filename);
  toTiffImage (noatoms, filename2);

  strstream3 << para << "_" << buffer << "_atom.raw";
  strstream3 >> filename;
  strstream4 << para << "_" << buffer << "_noatom.raw";
  strstream4 >> filename2;

  save_gsl_matrix_ASCII (filename, atoms);
  save_gsl_matrix_ASCII (filename2, noatoms);

  gsl_matrix_free (atoms);
  gsl_matrix_free (atomsNonCrop);

}
