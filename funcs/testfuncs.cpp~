/*
 * Project:  Testing Fermi-Dirac functions
 *
 * Author:   Pedro M Duarte 2011-02
 * 
 */

#include "funcs.h"
#include "math.h"

using namespace std;

bool VERBOSE;
bool DEBUG_FUNCS;

int
main (int argc, char **argv)
{

  VERBOSE = false;

  printf ("Fermi-Dirac  F0(-40.) = %.6e\n", f1 (-40.));
  printf ("Fermi-Dirac  F1(-40.) = %.6e\n", f1 (-40.));
  printf ("Fermi-Dirac  F2(-40.) = %.6e\n", f2 (-40.));
  printf ("Fermi-Dirac F32(-40.) = %.6e\n", f32 (-40.));
  printf ("Fermi-Dirac FM1(-40.) = %.6e\n", fm1 (-40.));
  printf ("Fermi-Dirac  FQ(-40.) = %.6e\n", fq (-40.));

  printf ("Fermi-Dirac  F1(-700) = %.6f\n", f1 (-700.));
  printf ("Fermi-Dirac  F1(-712.188276) = %.6f\n", f1 (-712.188276));
  /*for (int i = -100; i < 100; i++)
     printf ("%f\t%f\t%f\n", i * 1.0, f0 (i * 1.0) / fm1 (i * 1.0),
     fq (i * 1.0)); */
  //printf ("Fermi-Dirac  F1(-100.) = %.6f\n", f1 (-100.));
  return 0;
}
