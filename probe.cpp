/*
 * Project:  App3 Data analysis 
 *
 * File:     analyze.cpp
 *
 * Author:   Pedro M Duarte 2010-11
 * 
 */

#include "utils/utils.h"
#include "fits/fits.h"
#include "funcs/funcs.h"
#include <getopt.h>
#include <sstream>

#include "Fermions.h"

bool VERBOSE;


int processArgsAnalyze (int argc, char **argv, struct params &p);

int
main (int argc, char **argv)
{
  struct params p;
  init_params (&p);
  processArgsAnalyze (argc, argv, p);
  VERBOSE = p.verbose;

  Fermions *f = new Fermions (&p);
  f->LoadFITS ();
  f->ComputeColumnDensity ();


  printf
    ("%s  OD_max = %.2f, I_max = %.2f \n", p.shotnum.c_str (), f->maxOD,
     f->maxI);

  return EXIT_SUCCESS;


}


int
processArgsAnalyze (int argc, char **argv, struct params &p)
{
/*  Read command line arguments */


  if (argc == 1)
    {
      cout << endl;
      cout << "  usage:  " << argv[0] << " SHOTNUM [OPTIONS] " << endl;
      cout << endl;
      cout << "  OPTIONS " << endl << endl;

      cout << "\t-C, --center [c0,c1]"
	<< endl <<
	"\t\tmanuallly specify the initial guess for the cloud center" << endl
	<< "\t\t(not implemented yet, doesn't do anything)" << endl << endl;

      cout << "\t-c, --crop" << endl <<
	"\t\tcrop images before any further analysis" << endl << endl;

      cout << "\t-f, --force" << endl
	<< "\t\tforce reanalysis of the shot" << endl << endl;

      cout << "\t-p, --plots" << endl
	<< "\t\tplots profiles when fitting" << endl << endl;

      cout << "\t-r [PATH], --ref [PATH]" <<
	endl << "\t\tindicates path of reference image" << endl << endl;

      cout << "\t-R, --roi [ax0pos,ax1pos,ax0size,ax1size]"
	<< endl << "\t\tsets the atoms region of interest" << endl << endl;

      cout << "\t-S, --roisize [ax0size,ax1size]"
	<< endl <<
	"\t\tsets the size for the atoms region of interest" << endl <<
	"\t\tprogram attempts to center ROI around cloud" << endl << endl;

      cout << "\t-v, --verbose" << endl
	<< "\t\tshow messages to explain what is being done" << endl << endl;

      exit (2);
    }

  makeShotPaths (argv[1], p.shotnum, p.reportfile, p.atomsfile,
		 p.noatomsfile, p.atomsreffile, p.noatomsreffile);

  p.shot = atoi (p.shotnum.c_str ());
  p.verbose = false;
  p.reanalyze = false;
  p.center = false;
  p.crop = false;
  p.plots = false;
  p.roi_user = false;
  p.roisize_user = false;

  int c;
  string temp;
  stringstream ss (stringstream::in | stringstream::out);

  while (1)
    {
      static struct option long_options[] = {
	{"center", required_argument, 0, 'C'},
	{"crop", no_argument, 0, 'c'},
	{"force", no_argument, 0, 'f'},
	{"ref", required_argument, 0, 'r'},
	{"plots", no_argument, 0, 'p'},
	{"roi", required_argument, 0, 'R'},
	{"roisize", required_argument, 0, 'S'},
	{"verbose", no_argument, 0, 'v'},
	{0, 0, 0, 0}
      };


      int option_index = 0;
      c =
	getopt_long (argc, argv, "Ccfpr:R:S:v", long_options, &option_index);
      if (c == -1)
	break;

      switch (c)
	{
	case 'C':
	  p.center = 1;
	  temp = optarg;
	  replace (temp.begin (), temp.end (), ',', '\n');
	  ss << temp;
	  ss >> (p.centerpt)[0];
	  ss >> (p.centerpt)[1];
	  break;

	case 'c':
	  p.crop = 1;
	  break;

	case 'f':
	  p.reanalyze = 1;
	  break;

	case 'p':
	  p.plots = true;
	  break;

	case 'r':
	  p.atomsreffile = optarg;
	  p.noatomsreffile = optarg;
	  break;

	case 'R':
	  p.roi_user = true;
	  temp = optarg;
	  replace (temp.begin (), temp.end (), ',', '\n');
	  ss << temp;
	  ss >> (p.roi)[0];
	  ss >> (p.roi)[1];
	  ss >> (p.roi)[2];
	  ss >> (p.roi)[3];
	  break;

	case 'S':
	  p.roisize_user = true;
	  temp = optarg;
	  replace (temp.begin (), temp.end (), ',', '\n');
	  ss << temp;
	  ss >> (p.roisize)[0];
	  ss >> (p.roisize)[1];
	  break;

	case 'v':
	  p.verbose = 1;
	  break;

	case '?':
	  break;

	default:
	  abort ();
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
//  atomsfile << endl << p.noatomsfile << endl;

  /************* PARAMETERS OBTAINED FROM REPORT *******/
  p.texp = (double) getINI_num (p.reportfile, "ANDOR", "exp") * 1e-3;	// Exposure in seconds


  return EXIT_SUCCESS;
}
