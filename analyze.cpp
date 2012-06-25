/*
 * Project:  App3 Data analysis 
 *
 * File:     analyze.cpp
 *
 * Author:   Pedro M Duarte 2010-11
 * 
 */

#include "utils/utils.h"
#include "funcs/funcs.h"
#include "fits/fits.h"

#include "vt100_macros.h"
#include <getopt.h>
#include <time.h>
#include <sstream>

#include "Fermions.h"

bool VERBOSE;


int processArgsAnalyze (int argc, char **argv, struct params &p);
int writelog (int argc, char **argv);

int
main (int argc, char **argv)
{
  struct params p;
  processArgsAnalyze (argc, argv, p);
  init_params (&p);
  VERBOSE = p.verbose;


  Fermions *f = new Fermions (&p);
  f->LoadFITS ();		// LoadFITS already computes the column density

  setINI_num (p.reportfile, "CPP", "nsp", f->Nsp);
  setINI_num (p.reportfile, "CPP", "maxOD", f->maxOD);
  setINI_num (p.reportfile, "CPP", "maxPHI", f->maxPHI);
  setINI_num (p.reportfile, "CPP", "maxCD", f->maxCD);
  setINI_num (p.reportfile, "CPP", "maxI", f->maxI);

  f->FindMoments ();
  if (!p.keeproi)
    {
      f->MinimalCrop (5.0);
      f->FindMoments ();
    }

  f->Fit2DGauss ();
  f->SaveColumnDensity ();
  if (!p.keeproi)
    {
      f->MinimalCrop (3.5);
      f->FindMoments ();
      f->Fit2DGauss ();
      f->SaveColumnDensity ();
    }
  f->FitScatt2DGauss ();
  //f->FitProbe2DGauss ();
  f->Get2DCuts (true, false);

  setINI_num (p.reportfile, "CPP", "offset",
	      f->gaus2dfit[5] * f->GetNPixels ());

  setINI_num (p.reportfile, "CPP", "nfit", f->nfit);
  setINI_num (p.reportfile, "CPP", "nfit_err", f->nfit_err);
  setINI_num (p.reportfile, "CPP", "peakd", f->peakd);
  setINI_num (p.reportfile, "CPP", "ax0w", f->wi_1e * p.magnif);
  setINI_num (p.reportfile, "CPP", "ax1w", f->wj_1e * p.magnif);
  setINI_num (p.reportfile, "CPP", "ax0c", f->abs_ci);
  setINI_num (p.reportfile, "CPP", "ax1c", f->abs_cj);
  setINI_num (p.reportfile, "CPP", "TF", f->TF);
  setINI_num (p.reportfile, "CPP", "ph_per_at", f->Tsp / f->nfit);
  setINI_num (p.reportfile, "CPP", "peak_cd", f->peak);
  setINI_num (p.reportfile, "CPP", "photons_in_pulse", f->Tp0);
  setINI_num (p.reportfile, "CPP", "alphastar", p.alphastar);


  if (p.fermi2d)
    {
      f->Fit2DFermi ();
      setINI_num (p.reportfile, "CPP", "n0_Fermi", f->n0);
      setINI_num (p.reportfile, "CPP", "BetaMU_Fermi", f->BetaMu);
      setINI_num (p.reportfile, "CPP", "ri_Fermi", f->ri_Fermi);
      setINI_num (p.reportfile, "CPP", "rj_Fermi", f->rj_Fermi);
      setINI_num (p.reportfile, "CPP", "ci_Fermi", f->ci_Fermi);
      setINI_num (p.reportfile, "CPP", "cj_Fermi", f->cj_Fermi);
      setINI_num (p.reportfile, "CPP", "B_Fermi", f->B_Fermi);
      setINI_num (p.reportfile, "CPP", "TF_2d", f->TF_2d);
      setINI_num (p.reportfile, "CPP", "T_2d_rd", f->T_2d_rd);
      setINI_num (p.reportfile, "CPP", "T_2d_ax", f->T_2d_ax);
      f->Get2DCuts (true, p.fermi2d);
    }


  f->GetAzimuthalAverageEllipse ();

  if (p.fermiazimuth)
    {
      //f->GetAzimuthalAverage ();
      f->FitAzimuthalFermi ();
      setINI_num (p.reportfile, "CPP", "n0_az", f->n0_az);
      setINI_num (p.reportfile, "CPP", "BetaMu_az", f->BetaMu_az);
      setINI_num (p.reportfile, "CPP", "r_az", f->r_az);
      setINI_num (p.reportfile, "CPP", "B_az", f->B_az);
      setINI_num (p.reportfile, "CPP", "mx_az", f->mx_az);
      setINI_num (p.reportfile, "CPP", "TF_az", f->TF_az);
      setINI_num (p.reportfile, "CPP", "T_az", f->T_az);
    }

  if (p.fitfermi1D || p.fermiazimuth || p.fermi2d || true)
    {
      f->ComputeIntegrated1DDensity ();
      f->MakePlots ();
    }

  printf ("%s  N = %.2e", p.shotnum.c_str (), f->nfit);

  //Print number
  if (f->nfit_err * 1e2 > f->nfit)
    {
      printf
	("\nNumber determination uncertainty might be too high: N = %.2e +/- %.0e\n",
	 f->nfit, f->nfit_err);
    }

  //Get centr of cloud with respect to the Andor full frame
  f->abs_ci += f->gaus2dfit[0];
  f->abs_cj += f->gaus2dfit[2];


  if (!p.phc)
    {
      //Print number from scattered photons
      printf (", Nsp = %.2e", f->Nsp);

      //Print total number of photons in probe pulse
      printf (", Ph = %.2e", f->Tp0);

      //Print average number of photons scattered per atom
      printf (", Ph/At = %.0f", f->Tsp / f->nfit);

      //Print peak density
      printf (", n = %.2e", f->peakd);

      //Print axial size
      printf (", ax0w = %.1f +/- %.1f", f->wi_1e * p.magnif,
	      f->gaus2dfit_err[1] * p.magnif);

      //Print center
      printf (", c = (%.0f,%.0f)", f->abs_ci, f->abs_cj);

      //Print max intensity
      printf (", Imax = %.2f", f->maxI);

      //Print max optical density
      printf (", ODmax = %.1f", f->maxOD);

      //Print max column density
      printf (", CDmax = %.1f", f->maxCD);

    }

  else
    {
      //Print total number of photons in probe pulse
      printf (", Ph = %.2e", f->Tp0);

      //Print peak density
      printf (", n = %.2e", f->peakd);

      //Print axial size
      printf (", ax0w = %.1f +/- %.1f", f->wi_1e * p.magnif,
	      f->gaus2dfit_err[1] * p.magnif);

      //Print center
      printf (", c = (%.0f,%.0f)", f->abs_ci, f->abs_cj);

      //Print max intensity
      printf (", Imax = %.2f", f->maxI);

      //Print max phase shift
      printf (", PHImax = %.3f", f->maxPHI);

      //Print max column density
      printf (", CDmax = %.1f", f->maxCD);

    }
//  printf
//    ("%s  N = %.2e +/- %.0e, n =  %.2e , ax0w = %.1f +/- %.1f, c = (%.0f,%.0f), I_max = %.2f, OD_max = %.1f",
//     p.shotnum.c_str (), f->nfit, f->nfit_err, f->peakd,
//     f->wi_1e * p.magnif, f->gaus2dfit_err[1] * p.magnif, f->abs_ci,
//     f->abs_cj, f->maxI, f->maxOD);

  if (p.fermi2d)
    printf (", T/TF_2d = %.2f, z = %.2f, f(z) = %.2f", f->TF_2d,
	    f->Fugacity_Fermi, f->f_Fermi);
  if (p.fermiazimuth)
    printf (", T/TF_az = %.2f, z = %.2f, f(z) = %.2f", f->TF_az,
	    f->Fugacity_az, f->f_az);
  printf ("\n");


/*

  double pos[2];

  double cts = img_counts (signal);
  double peak = img_peak (signal, pos);

  setINI_num (p.reportfile, "CPP", "peak", peak);
  setINI_num (p.reportfile, "CPP", "ipeak", pos[0]);
  setINI_num (p.reportfile, "CPP", "jpeak", pos[1]);
  cout << "#" << p.
    shotnum << " Counts=" << img_counts (signal) << " Peak=" << peak << endl;
*/
  return EXIT_SUCCESS;

}


int
processArgsAnalyze (int argc, char **argv, struct params &p)
{
/*  Read command line arguments */
  writelog (argc, argv);

  if (argc == 1)
    {
      printf ("\n  usage:  %s SHOTNUM [OPTIONS]\n\n", argv[0]);
      printf (BOLDWHITE "  OPTIONS\n\n");

      printf (BOLDWHITE "\t--phc\n" RESET);
      printf ("\t\tdo polarization phase contrast analysis\n\n");

      printf (BOLDWHITE "\t--show-B\n" RESET);
      printf ("\t\tshow odttof, B1(w,t), and B2(w,t) then exit\n\n");

      printf (BOLDWHITE "\t-a, --alphastar [a*]\n" RESET);
      printf ("\t\tset absorption imaging calibration parameter\n\n");

      printf (BOLDWHITE "\t--trapfreq [vr]\n" RESET);
      printf ("\t\tuse to override default radial trapfreq in kHz \n\n");

      printf (BOLDWHITE "\t-C, --center [c0,c1]\n" RESET);
      printf
	("\t\tmanually specify the initial guess for the cloud center\n");
      printf ("\t\t(not implemented yet, does not do anything\n\n");

      printf (BOLDWHITE "\t-c, --crop\n" RESET);
      printf
	("\t\tcrop images, according to user speficied region before doing any fits\n\n");
      printf
	("\t\tNOTE: by default images are autocropped after a first fit with a 2D Gaussian.\n\n");
      printf
	("\t\t      If you use to keep the user defined ROI use --keeproi\n\n");

      printf (BOLDWHITE "\t--keeproi\n" RESET);
      printf ("\t\tkeeps the user defined ROI, does not autocrop\n\n");

      printf (BOLDWHITE "\t--fermi1d\n" RESET);
      printf ("\t\tperform fermi fits on integrated 1D profiles\n\n");

      printf (BOLDWHITE "\t--fermi2d\n" RESET);
      printf ("\t\tperform 2D Fermi-Dirac fit\n\n");

      printf (BOLDWHITE "\t--fermi-azimuth\n" RESET);
      printf ("\t\tperform azimuthal Fermi-Dirac fit\n\n");

      printf (BOLDWHITE "\t--maxr-azimuth [MAXR]\n" RESET);
      printf ("\t\tset max radius to be considered in azimuthal fits\n");
      printf ("\t\tif this option is given --chop-azimuth will be ignored\n");

      printf (BOLDWHITE "\t--chop-azimuth [DIST]\n" RESET);
      printf
	("\t\tdistance to be chopped off the tail when doing an azimuthal average\n");
      printf
	("\t\tthe tail is always noise as there are less points to do averaging with\n\n");

      printf (BOLDWHITE "\t--start-azimuth [DIST]\n" RESET);
      printf
	("\t\tdistance from the center for the start of  the azimuthal fit\n");
      printf
	("\t\tthe center is unwanted because the Fermi character of the cloud is more pronounced on the wings\n\n");

      printf (BOLDWHITE "\t--show-fermi\n" RESET);
      printf ("\t\tshow results of 2D Fermi and/or azimuthal Fermi fits\n\n");

      printf (BOLDWHITE "\t-f, --force\n" RESET);
      printf ("\t\tforce reanalysis of the shot\n\n");

      printf (BOLDWHITE "\t-p, --plots\n" RESET);
      printf ("\t\tplots profiles when fitting\n\n");

      printf (BOLDWHITE "\t-r [PATH], --ref [PATH]\n" RESET);
      printf ("\t\tindicates path of reference image\n\n");

      printf (BOLDWHITE "\t-R, --roi [ax0pos,ax1pos,ax0size,ax1size]\n"
	      RESET);
      printf ("\t\tsets the atoms region of interest\n\n");

      printf (BOLDWHITE "\t-S, --roisize [ax0size,ax1size]\n" RESET);
      printf ("\t\tsets the size for the atoms region of interest\n");
      printf
	("\t\tthe program will automatically determine the center for this ROI\n\n");

      printf (BOLDWHITE "\t--blanks\n" RESET);
      printf
	("\t\tuse this option when taking empty pictures (for diagnosing probe, etc.)\n\n");

      printf (BOLDWHITE "\t-v, --verbose\n" RESET);
      printf ("\t\tshow messages to explain what is being done\n\n");

      exit (2);
    }

  makeShotPaths (argv[1], p.shotnum, p.reportfile, p.atomsfile,
		 p.noatomsfile, p.atomsreffile, p.noatomsreffile);

  p.shot = atoi (p.shotnum.c_str ());
  p.verbose = false;
  p.reanalyze = false;
  p.center = false;
  p.crop = false;
  p.keeproi = false;
  p.plots = false;
  p.roi_user = false;
  p.roisize_user = false;
  p.fitfermi1D = false;
  p.fermi2d = false;
  p.fermiazimuth = false;
  p.azimuth_maxr = 512.;
  p.azimuth_chop = 0.;
  p.azimuth_start = 0.;
  p.showfermi = false;
  p.show_B = false;		// Use this to show B factors, which affect temperature determination vio cloud size
  p.w_user = false;		// this flag is set if the user specifies the radial trap frequency from the command line 

  p.blanks = false;		// this flag us used when the pictures are blanks (i.e. have no atoms) it is useful for probe diagnostics
  // notice that this flag can be set on the command line or captured from the report.

  p.alphastar = 1.0;		// Calibration for absorption imaging polarizaation and effective saturation
  // See "Strong saturation of absorption imaging of dense clouds of ultracold atoms" 
  //      G. Reinaudi et al.   Optics Letters, Vol. 32, No. 21 pg 3143  
  p.phc = false;

  int c;
  string temp;
  stringstream ss (stringstream::in | stringstream::out);

  while (1)
    {
      static struct option long_options[] = {
	{"alphastar", required_argument, 0, '*'},
	{"show-B", no_argument, 0, 'B'},
	{"phc", no_argument, 0, 'P'},
	{"center", required_argument, 0, 'C'},
	{"crop", no_argument, 0, 'c'},
	{"keeproi", no_argument, 0, 'k'},
	{"fermi1d", no_argument, 0, '+'},
	{"fermi2d", no_argument, 0, 'F'},
	{"fermi-azimuth", no_argument, 0, 'a'},
	{"maxr-azimuth", required_argument, 0, 'd'},
	{"chop-azimuth", required_argument, 0, 'h'},
	{"start-azimuth", required_argument, 0, 't'},
	{"show-fermi", no_argument, 0, 's'},
	{"force", no_argument, 0, 'f'},
	{"ref", required_argument, 0, 'r'},
	{"plots", no_argument, 0, 'p'},
	{"roi", required_argument, 0, 'R'},
	{"roisize", required_argument, 0, 'S'},
	{"blanks", no_argument, 0, 'b'},
	{"verbose", no_argument, 0, 'v'},
	{"trapfreq", required_argument, 0, 'w'},
	{0, 0, 0, 0}
      };


      int option_index = 0;
      c =
	getopt_long (argc, argv, "Ccfpr:R:S:v", long_options, &option_index);
      if (c == -1)
	break;

      switch (c)
	{
	case '*':
	  temp = optarg;
	  ss << temp;
	  ss >> p.alphastar;

	  break;

	case 'P':
	  p.phc = true;
	  break;

	case 'B':
	  p.show_B = true;
	  break;

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

	case 'k':
	  p.keeproi = true;
	  break;

	case '+':
	  p.fitfermi1D = true;
	  break;

	case 'F':
	  p.fermi2d = true;
	  break;

	case 'a':
	  p.fermiazimuth = true;
	  break;

	case 'd':
	  p.azimuth_maxr = atof (optarg);
	  break;

	case 'h':
	  p.azimuth_chop = atof (optarg);
	  break;

	case 't':
	  p.azimuth_start = atof (optarg);
	  break;

	case 's':
	  p.showfermi = true;
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

	case 'b':
	  p.blanks = true;
	  break;

	case 'v':
	  p.verbose = 1;
	  break;

	case 'w':
	  p.w_user = true;
	  temp = optarg;
	  ss << temp;
	  ss >> p.w;
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
  p.odttof = (double) getINI_num (p.reportfile, "ODT", "odttof");	// odttof in ms 
  p.det = (double) getINI_num (p.reportfile, "ANDOR", "phcdet");	// phase contrast detuning in MHz
  p.phc = (bool) getINI_num (p.reportfile, "ANDOR", "phc");	// phase contrast 
  p.phcsign = (double) getINI_num (p.reportfile, "ANDOR", "phcsign");	// prefactor in the calculation of the phase contrast signal

  bool auxblanks = (bool) getINI_num (p.reportfile, "ANDOR", "blanks");	// taking blanks ?
  p.blanks = p.blanks || auxblanks;	// if either the command line or the report says it, these are blanks

//  p.det = -240.;
//  p.phc = 0.;
//  p.phcsign = 1.;
  double magnif_INI = (double) getINI_num (p.reportfile, "ANDOR", "magnif");
  if (magnif_INI == INI_ERROR)
    printf
      ("\n   Did not find ANDOR:magnif in the report file.\n   Using hard-coded value %.2f um/pixel\n\n",
       p.magnif);
  else
    p.magnif = magnif_INI;

  p.free = (double) getINI_num (p.reportfile, "EVAP", "free");	// 
  p.p1 = (double) getINI_num (p.reportfile, "EVAP", "p1");	// 
  p.t1 = (double) getINI_num (p.reportfile, "EVAP", "t1");	// 
  p.tau = (double) getINI_num (p.reportfile, "EVAP", "tau");	// 
  p.beta = (double) getINI_num (p.reportfile, "EVAP", "beta");	// 
  p.p0 = (double) getINI_num (p.reportfile, "ODT", "odtpow0");	// 
  p.offset = (double) getINI_num (p.reportfile, "EVAP", "offset");	// 
  p.t2 = (double) getINI_num (p.reportfile, "EVAP", "t2");	// 
  p.tau2 = (double) getINI_num (p.reportfile, "EVAP", "tau2");	// 
  p.image = (double) getINI_num (p.reportfile, "EVAP", "image");	// 


  return EXIT_SUCCESS;
}


int
writelog (int argc, char **argv)
{

  time_t rawtime;
  struct tm *timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime (&rawtime);

  strftime (buffer, 80, " %a %d %b %Y %X %p %Z   |   ", timeinfo);



  ofstream fout ("analysislog", ofstream::app);
//   fout << "analyze  " ;
  fout << buffer;
  for (int i = 0; i < argc; i++)
    {
      fout << argv[i] << " ";
    }
  fout << endl;
  fout.close ();
  return 0;
}
