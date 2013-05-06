/*
 * Project:  App3 Data analysis 
 *
 * File:     analyze.cpp
 *
 * Author:   Pedro M Duarte 2010-11
 * 
 */

//Some shorthand to output colors to terminal
#include <stdio.h>
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"


#include "utils/utils.h"
#include "funcs/funcs.h"
#include "fits/fits.h"

#include "vt100_macros.h"
#include <getopt.h>
#include <time.h>
#include <sstream>

#include "Fermions.h"





bool VERBOSE;
bool DEBUG_FUNCS;
bool DEBUG_FITS;


int processArgsAnalyze (int argc, char **argv, struct params &p);
int writelog (int argc, char **argv);

int
main (int argc, char **argv)
{
  struct params p;
  processArgsAnalyze (argc, argv, p);
  VERBOSE = p.verbose;
  init_params (&p);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "phcdet", p.det);


  Fermions *f = new Fermions (&p);
  f->LoadFITS ();		// LoadFITS already computes the column density

  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "maxOD", f->maxOD);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "maxPHI", f->maxPHI);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "maxPHCSIG",
	      f->maxPHCSIG);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "maxCD", f->maxCD);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "maxI", f->maxI);

  if (p.onlyCD)
    {
      f->SaveColumnDensity ();
      return EXIT_SUCCESS;
    }
//  f->FindMoments ();
  if (!p.keeproi)
    {
      f->MinimalCrop (5.0);
    }
  f->Fit2DGauss ();
  f->SaveColumnDensity ();
  if (!p.keeproi)
    {
      f->MinimalCrop (3.5);
      f->Fit2DGauss ();
      f->SaveColumnDensity ();
    }
//  f->FitScatt2DGauss ();
//  f->FitProbe2DGauss ();

  //Get center of cloud with respect to the Andor full frame
  f->abs_ci += f->gaus2dfit[0];
  f->abs_cj += f->gaus2dfit[2];

  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "offset",
	      f->gaus2dfit[5] * f->GetNPixels ());

  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "nfit", f->nfit);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "nfit_err",
	      f->nfit_err);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "peakd", f->peakd);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "peakd_sph",
	      f->peakd_sph);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "ax0w",
	      f->gaus2dfit[1] * p.magnif);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "ax1w",
	      f->gaus2dfit[3] * p.magnif);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "ax0c", f->abs_ci);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "ax1c", f->abs_cj);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "TF", f->TF);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "ph_per_at",
	      f->Tsp / f->nfit);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "peak_cd",
	      f->gaus2dfit[4]);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "photons_in_pulse",
	      f->Tp0);
  setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "alphastar",
	      p.alphastar);

  f->GetAzimuthalAverageEllipse ();

  if (p.fermiazimuth)
    {
      //f->GetAzimuthalAverage ();
      f->FitAzimuthalFermi ();
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "n0_az", f->n0_az);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "BetaMu_az",
		  f->BetaMu_az);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "r_az", f->r_az);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "B_az", f->B_az);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "mx_az", f->mx_az);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "TF_az", f->TF_az);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "T_az", f->T_az);
    }

  if (p.fermi2d)
    {
      f->Fit2DFermi ();
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "n0_Fermi", f->n0);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "BetaMU_Fermi",
		  f->BetaMu);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "ri_Fermi",
		  f->ri_Fermi);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "rj_Fermi",
		  f->rj_Fermi);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "ci_Fermi",
		  f->ci_Fermi);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "cj_Fermi",
		  f->cj_Fermi);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "B_Fermi",
		  f->B_Fermi);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "TF_2d", f->TF_2d);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "T_2d_rd",
		  f->T_2d_rd);
      setINI_num (p.reportfile, p.andor2 ? "CPP2" : "CPP", "T_2d_ax",
		  f->T_2d_ax);
    }


  if (p.fitfermi1D || p.fermiazimuth || p.fermi2d || true)
    {
      f->ComputeIntegrated1DDensity ();
      f->MakePlots ();
    }

  printf ("%s  %sN = %.2e%s", p.shotnum.c_str (), KCYN, f->nfit, KNRM);

  //Print number
  if (f->nfit_err * 50. > f->nfit)
    {
      printf
	("\nNumber determination uncertainty might be too high: N = %.2e +/- %.0e\n",
	 f->nfit, f->nfit_err);
    }



  if (!p.phc && !p.fluor)
    {
      //Print number from scattered photons
      //printf (", Nsp = %.2e", f->Nsp);

      //Print total number of photons in probe pulse
      //printf (", Ph = %.2e", f->Tp0);

      //Print average number of photons scattered per atom
      printf (", Ph/At = %.0f", f->Tsp / f->nfit);

      //Print peak density
      printf (", n = %.2e", f->peakd);

      //Print peak density
      printf (", %sn_sph = %.2e%s", KMAG, f->peakd_sph, KNRM);

      //Print axial size
      printf (", ax0w = %.1f +/- %.1f", f->gaus2dfit[1] * p.magnif,
	      f->gaus2dfit_err[1] * p.magnif);

      //Print center
      printf (", %sc = (%.0f,%.0f)%s", KGRN, f->abs_ci, f->abs_cj, KNRM);

      //Print max intensity
      //printf (", Imax = %.2f", f->maxI);

      //Print max intensity smoothed
      //printf (", Imax(smoothed) = %.2f", f->maxIsmooth);

      //Print average intensity
      //printf (", Iave = %.2f", f->aveI);

      //Print average intensity
      printf (", I/Isat = %.2f", f->aveIweighted);

      //Print max optical density
      printf (", ODmax = %.1f", f->maxOD);
      //Print max column density
      printf (", CDmax = %.1f", f->maxCD);

    }

  else if (p.fluor || p.Nframes == 2)
    {
      printf (", MAX Counts/Px = %.1f", f->maxCD);

      //Print axial size
      printf (", ax0w = %.1f +/- %.1f", f->gaus2dfit[1] * p.magnif,
	      f->gaus2dfit_err[1] * p.magnif);

      //Print center
      printf (", c = (%.0f,%.0f)", f->abs_ci, f->abs_cj);

      //Print range of counts in atoms and noatoms frames
      printf (", atoms:(%d to %d),  noatoms:(%d to %d)", (int) f->minCA,
	      (int) f->maxCA, (int) f->minCN, (int) f->maxCN);
    }

  else
    {
      //Print total number of photons in probe pulse
      printf (", Ph = %.2e", f->Tp0);
      //Print peak density
      printf (", n = %.2e", f->peakd);
      //Print peak density
      printf (", %sn_sph = %.2e%s", KMAG, f->peakd_sph, KNRM);
      //Print axial size
      printf (", ax0w = %.1f +/- %.1f", f->gaus2dfit[1] * p.magnif,
	      f->gaus2dfit_err[1] * p.magnif);
      //Print center
      printf (", %sc = (%.0f,%.0f)%s", KGRN, f->abs_ci, f->abs_cj, KNRM);
      //Print max intensity
      //printf (", Imax = %.2f", f->maxI);
      //Print average intensity
      printf (", I/Isat = %.2f", f->aveIweighted);
      //Print max phase shift
      printf (", PHImax = %.3f", f->maxPHI);
      //Print max phase-contrast signal
      printf (", SIGmax = %.3f", f->maxPHCSIG);
      //Print max column density
      printf (", CDmax = %.1f", f->maxCD);

      //Print sqrt(N)/n parameter
      printf (", %ssqrt(N)/n = %.2e%s", KYEL, sqrt (f->nfit) / f->peakd_sph,
	      KNRM);

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

  setINI_num (p.reportfile, p.andor2?"CPP2":"CPP", "peak", peak);
  setINI_num (p.reportfile, p.andor2?"CPP2":"CPP", "ipeak", pos[0]);
  setINI_num (p.reportfile, p.andor2?"CPP2":"CPP", "jpeak", pos[1]);
  cout << "#" << p.
    shotnum << " Counts=" << img_counts (signal) << " Peak=" << peak << endl;
*/
  return EXIT_SUCCESS;
}


int
processArgsAnalyze (int argc, char **argv, struct params &p)
{
/*  Read command line arguments */
  char argv1[MAXPATHLEN];
  if (argc > 1)
    strcpy (argv1, argv[1]);

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

      printf (BOLDWHITE "\t--cdsp\n" RESET);
      printf
	("\t\tuses number from scattered photons to obtain the column density\n\n");

      printf (BOLDWHITE "\t--highintabs [PROBEPOWER]\n" RESET);
      printf
	("\t\tuses high intensity number from scattered photons to obtain the column density\n");
      printf ("\t\tthe power of the probe beam must be provided. \n\n");

      printf (BOLDWHITE "\t--magnif [m]\n" RESET);
      printf ("\t\tuse to override the magnification. m is in um/pixel\n\n");

      printf (BOLDWHITE "\t--probewaist [w]\n" RESET);
      printf ("\t\tuse to override the probe beam waist. w is in cm\n\n");

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
      printf
	("\t\tif this option is given --chop-azimuth will be ignored\n\n");

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

      printf (BOLDWHITE "\t-r [PATH], --ref [PATH]\n" RESET);
      printf ("\t\tindicates path of reference image\n\n");

      printf (BOLDWHITE
	      "\t-R, --roi [ax0pos,ax1pos,ax0size,ax1size]\n" RESET);
      printf ("\t\tsets the atoms region of interest\n\n");

      printf (BOLDWHITE "\t--blanks\n" RESET);
      printf
	("\t\tuse this option when taking empty pictures (for diagnosing probe, etc.)\n\n");

      printf (BOLDWHITE "\t--saveascii\n" RESET);
      printf
	("\t\tuse this option to enable saving of the ascii files with image processing data.\n\n");

      printf (BOLDWHITE "\t--savetiff\n" RESET);
      printf
	("\t\tuse this option to enable saving of the TIFF files with image processing data.\n\n");

      printf (BOLDWHITE "\t-v, --verbose\n" RESET);
      printf ("\t\tshow messages to explain what is being done\n\n");

      printf (BOLDWHITE "\t-i, --imgverbose\n" RESET);
      printf
	("\t\tshow messages to explain the results of the column density calculation\n\n");

      printf (BOLDWHITE "\t--azimverbose\n" RESET);
      printf
	("\t\tshow histogram calculation for azimuthal averaging (produces long output)\n\n");

      printf (BOLDWHITE "\t-i, --onestate\n" RESET);
      printf
	("\t\tuse this flag if taking pictures of only state |1> atoms\n\n");

      printf (BOLDWHITE "\t-e, --eigenface\n" RESET);
      printf
	("\t\tuse this flag if you want to clean up the image with eigenface\n\n");

      printf (BOLDWHITE "\t--andor2\n" RESET);
      printf
	("\t\tuse this to analyze the set of pictures taken by andor2\n\n");

      printf (BOLDWHITE "\t--onlyCD\n" RESET);
      printf
	("\t\tuse this to only compute the column density, does not do any fits\n\n");

//      printf (BOLDWHITE "\t--debug-fits\n" RESET);
//      printf ("\t\tshow details about fit evaluations for debugging\n\n");
      exit (2);
    }

  p.shot = atoi (p.shotnum.c_str ());
  p.verbose = false;
  p.imgverbose = false;
  p.azimverbose = false;
  p.reanalyze = false;
  p.center = false;
  p.crop = false;
  p.keeproi = false;
  p.roi_user = false;
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
  p.twostates = true;		// By default pictures are of a spin mixture
  p.eigenface = false;		// By default do not use eigenface
  p.saveascii = false;
  p.savetiff = false;
  //Obtain column density using number from scattered photons in the high intensity limit
  p.highintabs = false;
  p.magnif_override = false;
  p.probewaist_override = false;
  //Obtain colulmn density using number from scattered photons
  p.cdsp = false;

  p.andor2 = false;		//By default andor2 analysis is false

  p.onlyCD = false;		//By default we want to fit the column density to something

  int c;
  while (1)
    {
      static struct option long_options[] = {
	{
	 "alphastar", required_argument, 0, '*'},
	{
	 "show-B", no_argument, 0, 'B'},
	{
	 "phc", no_argument, 0, 'P'},
	{
	 "center", required_argument, 0, 'C'},
	{
	 "crop", no_argument, 0, 'c'},
	{
	 "keeproi", no_argument, 0, 'k'},
	{
	 "fermi1d", no_argument, 0, '+'},
	{
	 "fermi2d", no_argument, 0, 'F'},
	{
	 "fermi-azimuth", no_argument, 0, 'a'},
	{
	 "maxr-azimuth", required_argument, 0, 'd'},
	{
	 "chop-azimuth", required_argument, 0, 'h'},
	{
	 "start-azimuth", required_argument, 0, 't'},
	{
	 "show-fermi", no_argument, 0, 's'},
	{
	 "force", no_argument, 0, 'f'},
	{
	 "ref", required_argument, 0, 'r'},
	{
	 "roi", required_argument, 0, 'R'},
	{
	 "blanks", no_argument, 0, 'b'},
	{
	 "saveascii", no_argument, 0, 'A'},
	{
	 "savetiff", no_argument, 0, 'T'},
	{
	 "verbose", no_argument, 0, 'v'},
	{
	 "imgverbose", no_argument, 0, 'i'},
	{
	 "azimverbose", no_argument, 0, 'z'},
	{
	 "onestate", no_argument, 0, '1'},
	{
	 "eigenface", no_argument, 0, 'e'},
	{
	 "trapfreq", required_argument, 0, 'w'},
	{
	 "highintabs", required_argument, 0, 'H'},
	{
	 "cdsp", no_argument, 0, 'n'},
	{
	 "magnif", required_argument, 0, 'm'},
	{
	 "probewaist", required_argument, 0, 'p'},
	{
	 "andor2", no_argument, 0, '2'},
	{
	 "onlyCD", no_argument, 0, 'D'},
	{
	 0, 0, 0, 0}
      };
      int option_index = 0;
      c =
	getopt_long (argc, argv, "Ccfpr:R:S:vie", long_options,
		     &option_index);
      if (c == -1)
	break;
      string temp;
      stringstream ss (stringstream::in | stringstream::out);
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
	case 'b':
	  p.blanks = true;
	  break;
	case 'A':
	  p.saveascii = true;
	  break;
	case 'T':
	  p.savetiff = true;
	  break;
	case 'v':
	  p.verbose = 1;
	  break;
	case 'i':
	  p.imgverbose = true;
	  break;
	case 'z':
	  p.azimverbose = true;
	  break;
	case 'e':
	  p.eigenface = true;
	  break;
	case '2':
	  p.andor2 = true;
	  break;
	case '1':
	  p.twostates = false;
	  break;
	case 'w':
	  p.w_user = true;
	  temp = optarg;
	  ss << temp;
	  ss >> p.w;
	  break;
	case 'm':
	  p.magnif_override = true;
	  temp = optarg;
	  ss << temp;
	  ss >> p.magnif_user;
	  break;
	case 'p':
	  p.probewaist_override = true;
	  temp = optarg;
	  ss << temp;
	  ss >> p.probewaist_user;
	  break;
	case 'H':
	  p.highintabs = true;
	  temp = optarg;
	  ss << temp;
	  ss >> p.probepower;
	  break;
	case 'n':
	  p.cdsp = true;
	  break;
	case 'D':
	  p.onlyCD = true;
	  break;
	case '?':
	  break;
	default:
	  abort ();
	}
    }

  makeShotPaths (argv1, p.shotnum, p.reportfile, p.atomsfile, p.noatomsfile,
		 p.atomsreffile, p.noatomsreffile, p.andor2);

  if (sectionExists (p.reportfile, p.andor2 ? "CPP2" : "CPP") && !p.reanalyze)
    {
      cout << endl;
      cout << " Shot " << p.shotnum << " has already been analyzed." << endl;
      cout << " Use option -f  to force analysis." << endl << endl;
      exit (2);
    }
//  cout << p.shotnum << endl << p.reportfile << endl << p.
//  atomsfile << endl << p.noatomsfile << endl;




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
