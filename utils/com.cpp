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


int
processArgsAnalyze (int argc, char **argv, struct params &p)
{
/*  Read command line arguments */


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

bool VERBOSE;

int
main (int argc, char **argv)
{

  VERBOSE = 1;

  struct params p;
  processArgsAnalyze (argc, argv, p);
  gsl_matrix *atomsNonCrop = ReadFitsImg_gsl_matrix (p.atomsfile);
  //gsl_matrix *noatoms = ReadFitsImg_gsl_matrix (p.noatomsfile);


  gsl_matrix *atoms = cropImage (p.reportfile, atomsNonCrop);
  //gsl_matrix *cnoatoms = cropImage (p.reportfile, noatoms);

  double com[2];
  COM (atoms, com);

  cout << "COM_X = " << com[0] << endl;
  cout << "COM_Y = " << com[1] << endl;

  setINI_num (p.reportfile, "CPP", "COM_X", com[0]);
  setINI_num (p.reportfile, "CPP", "COM_Y", com[1]);

  gsl_matrix_free (atomsNonCrop);
  gsl_matrix_free (atoms);

}
