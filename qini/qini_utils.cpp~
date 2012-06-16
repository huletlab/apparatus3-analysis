
#include "qini.h"


double
getINI_num (string & inifile, char *SECTION, char *KEY)
{

  //Loads INI file
  CSimpleIni reportINI (false, false, false);
  if (reportINI.LoadFile (inifile.c_str ()) < 0)
    {
      cout << "#Failed to load " << inifile << endl;
      exit (1);
    }
  const char *valuepointer;
  valuepointer = reportINI.GetValue (SECTION, KEY, NULL /*default */ );
  string value;
  if (valuepointer != NULL)
    {
      value = valuepointer;
    }
  else
    {
      cout << "#" << SECTION << ":" << KEY << " does not exist!" << endl;
      return INI_ERROR;
      //exit (EXIT_FAILURE);
    }

  double val = strtod (value.c_str (), NULL);
  if (isnan (val))
    {
      cout << "#" << " Error the value for " << SECTION << ":" << KEY <<
	" is NaN!" << endl;
      exit (EXIT_FAILURE);
    }
  return val;
}

int
setINI_num (string & inifile, char *SECTION, char *KEY, double val)
{
  //Loads INI file
  CSimpleIni reportINI (false, false, false);
  if (reportINI.LoadFile (inifile.c_str ()) < 0)
    {
      cout << "#Failed to load " << inifile << endl;
      exit (1);
    }
  char valstr[256];
  sprintf (valstr, "%.6e", val);
  reportINI.SetValue (SECTION, KEY, valstr, NULL);
  //  Save the data back to the file
  SI_Error rc = reportINI.SaveFile (inifile.c_str ());
  if (rc < 0)
    {
      cout << "Failure saving report file" << endl;
      exit (EXIT_FAILURE);
    }

  return EXIT_SUCCESS;

}

bool
sectionExists (string & inifile, char *SECTION)
{
  //Loads INI file
  CSimpleIni reportINI (false, false, false);
  if (reportINI.LoadFile (inifile.c_str ()) < 0)
    {
      cout << "#Failed to load " << inifile << endl;
      exit (1);
    }
  if (reportINI.GetSectionSize (SECTION) >= 0)
    return true;
  return false;
}
