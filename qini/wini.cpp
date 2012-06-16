/*
 * Project:  Set data fields in a .INI file
 *
 * Contents: This program sets a single value in a .INI file, given a section:key.
 *
 * Author:   Pedro M Duarte, Yihan Huang  2012-05-08
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
  string inifile;


  char copyargv[argc][256];
  for (int i = 0; i < argc; i++)
    {
      strcpy (copyargv[i], argv[i]);
    }


// Read command line arguments. getopt will parse the command line and
// return all options (-X) and their arguments. The third parameter of optarg
// is the list of valid options, a colon indicates the preceding option takes an 
// argument.
  int nopts = 0;

  opterr = 0;

  while ((c = getopt (argc, argv, "hf:")) != -1)
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
	case '?':
	  //printf ("found an unrecognized option\n");
	  break;
	default:
	  break;
	}
    }

  if (errflg || argc == 1)
    {
      cout << endl;
      cout << "  usage: " << argv[0] <<
	" [OPTIONS] [Shotnum] [section:key] [value]" << endl;
      cout << endl;
      cout <<
	"    Can be used with a shot number in the current working directory\n    or with an absolute path"
	<< endl;
      cout << endl;
      cout << "  Examples:" << endl << endl;
      cout << "   -Use with Shotnum:" << endl;
      cout << "    wini 5 section:key value" << endl << endl;
      cout << "   -Use with FilePath:" << endl;
      cout << "    wini -f /lab/data/***  section:key value" << endl << endl;
      cout << "  wini -h for help " << endl << endl;
      exit (2);
    }
  struct data d;
  vector < string > keys;

  if (argc > 1 && errflg == false)
    {

      //first argument is the shot number
      int shot;
      char shot_str[4];

/*      for (int i = 0; i < argc; i++)
	{
	  cout << "argument # " << i << " = " << argv[i] << " was " <<
	    copyargv[i] << endl;
	}*/
      //cout << "nopts = " << nopts << endl;

      //creates report file path in current working directory
      if (file == false)
	{
	  shot = atoi (copyargv[1 + nopts]);
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
      //cout << "inifile = " << inifile << endl;

      ofstream myfile;
      myfile.open (inifile.c_str (), ios_base::app);
      myfile.close ();

      string section;
      string key;
      string value;

      string str (copyargv[2 + nopts]);

      size_t found;
      found = str.find (":");
      section = str.substr (0, found);
      key = str.substr (found + 1);
      value = copyargv[3 + nopts];


      setINI (inifile, section.c_str (), key.c_str (), value.c_str ());
    }
  return 0;
}
