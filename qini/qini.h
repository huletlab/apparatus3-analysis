
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <unistd.h>
#include "simpleini/SimpleIni.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include <map>

#define INI_ERROR -1e11

using namespace std;

double getINI_num( string & inifile, const char *SECTION, const char *KEY);
int setINI_num (string & inifile, const char *SECTION, const char *KEY, double val);
int setINI (string & inifile, const char *SECTION, const char *KEY, const char *val);
bool sectionExists( string & inifile, const char *SECTION);
