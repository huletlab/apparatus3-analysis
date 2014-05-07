/*
 * Project:  Extracting data fields from a .INI file
 *
 * Contents: This program extracts a list of values given a section:key list from a
 *           given .INI file
 *
 * Author:   Pedro M Duarte 2010-09-22
 * 
 */


#include "qini.h"


struct data
{
  string shotnum;
  string reportfile;
};

int
main (int argc, char **argv)
{
  int c;
  bool errflg = false;
  bool file = false;
  bool raw = false;
  string inifile;
// Read command line arguments. getopt will parse the command line and
// return all options (-X) and their arguments. The third parameter of optarg
// is the list of valid options, a colon indicates the preceding option takes an 
// argument.
  unsigned int nopts = 0;
  while ((c = getopt (argc, argv, "hf:r")) != -1)
    {
      switch (c)
	{
	case 'h':
	  errflg = true;
	  break;
	case 'f':
	  file = true;
	  inifile = optarg;
	  nopts++;
	  break;
	case 'r':
	  raw = true;
	  nopts++;
	  break;
	default:
	  break;
	}
    }

  if (errflg || argc == 1)
    {
      cout << endl;
      cout << "  usage: " << argv[0] <<
	" [OPTIONS] [Shotnum] [list of section:key]" << endl;
      cout << endl;
      cout <<
	"    Can be used with a shot number in the current working directory\n    or with an absolute path"
	<< endl;
      cout << endl;
      cout << "  Examples:" << endl << endl;
      cout << "   -Use with Shotnum:" << endl;
      cout << "    qini 5 section1:key1 section2:key2 ..." << endl << endl;
      cout << "   -Use with FilePath:" << endl;
      cout << "    qini -f /lab/data/***  section1:key1 section2:key2 ..." <<
	endl << endl;
      cout <<
	"  Add the -r (raw) option to obtain the value string with no headers."
	<< endl << endl;
      cout << "  qini -h for help " << endl << endl;
      exit (2);
    }
  struct data d;
  vector < string > keys;

  if (argc > 1 && errflg == false)
    {

      //first argument is the shot number
      int shot;
      char shot_str[4];


      //creates report file path in current working directory
      if (file == false)
	{
	  shot = atoi (argv[1 + nopts]);
	  sprintf (shot_str, "%04d", shot);
	  d.shotnum = shot_str;
	  inifile = "";
	  char path[MAXPATHLEN];
	  getcwd (path, MAXPATHLEN);
	  string tmp3 ("");
	  inifile = tmp3;
	  inifile += path;
	  inifile += "/report";
	  inifile += shot_str;
	  inifile += ".INI";
	}
      //cout << argc << endl;

      //Loads INI file
      CSimpleIni reportINI (false, false, false);
      if (reportINI.LoadFile (inifile.c_str ()) < 0)
	{
	  cout << "#Failed to load " << inifile << endl;
	  exit (1);
	}

      //stores all section:key in a vector
      int start;
      start = 2 + nopts;
      for (int i = start; i < argc; i++)
	{
	  keys.push_back (argv[i]);
	}

      //outputs header line
      ostringstream head (ostringstream::out);
      ostringstream data (ostringstream::out);


      head << "#";
      data << " ";

      string value;
      for (vector < string >::iterator i = keys.begin (); i != keys.end ();
	   ++i)
	{
	  string section (i->substr (0, i->find (":")));
	  string key (i->substr (i->find (":") + 1));
	  const char *valuepointer;
	  valuepointer = reportINI.GetValue (section.c_str (), key.c_str (),
					     NULL /*default */ );
	  if (valuepointer != NULL)
	    {
	      value = valuepointer;
	    }
	  else
	    {
	      cout << "#" << section.c_str () << ":" << key.
		c_str () << " does not exist!" << endl;
	      exit (1);
	    }
	  int w = max (key.length (), value.length ());
	  head << setw (w + 3);
	  head << key;
	  data << setw (w + 3);
	  data << value;
	}
      if (raw)
	{
	  cout << data.str () << endl;
	}
      else
	{
	  cout << head.str () << endl << data.str () << endl;
	}
    }



  return 0;
}
