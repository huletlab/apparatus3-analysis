/*
 * Project:  This file defines various functions to handle file
 *           paths and files in general
 *           
 *
 * Author:   Pedro M Duarte 2011-01
 * 
 */


#include "utils/utils.h"

string
makepath (const char *base, string prefix, const char *identifier)
{
  string out ("");
  out += base;
  out += "/";
  out += prefix;
  out += identifier;
  return out;
}

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
