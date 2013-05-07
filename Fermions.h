/*
 *
 *  Evaluator class for absorption images of Fermions taken with Apparatus 3
 *
 *
 *  Author : Pedro M Duarte  Feb 2011
 *
 * 
 */


// This structure contains all parameters relevant for analysis

// Convention is  size1 : 0 : radial : i
//                size2 : 1 : axial  : j 


#include <math.h>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

extern bool VERBOSE;


struct params
{
  unsigned int shot;

  string shotnum, reportfile, atomsfile, noatomsfile, atomsreffile,
    noatomsreffile;

  string shotnum_fileout;

  // ROI is defined as (ax0pos, ax1pos, ax0size, ax1size)
  unsigned int roi[4], roisize[2], centerpt[2];
  bool keeproi;
  int Nframes;

  bool verbose, center, crop, plots, reanalyze, roi_user,
    fitfermi1D, fermi2d, fermiazimuth, showfermi, show_B, blanks, gaus2d_mott;

  // Use high intensity number from scattered photons in absorption imaging
  bool highintabs;
  double probepower;		//power of probe beam is required for this calculation
  double probewaist;		//

  // Use column density from scattered photons in absorption imaging
  bool cdsp;

  // Override the magnification 
  bool magnif_override;
  double magnif_user;

  // Override the probe beam waist
  bool probewaist_override;
  double probewaist_user;

  // One state or two states?
  bool twostates;

  // Be verbose about imaging
  bool imgverbose;

  // Be verbose about calculating azimuthal averages
  bool azimverbose;

  // Do eigenface cleanup?
  bool eigenface;
  bool eigenface_done;

  // Whether or not to save ASCII or TIFF files
  // Useful for debugging the algorithms
  bool saveascii, savetiff;


  double azimuth_maxr, azimuth_chop, azimuth_start;

  double lambda, hc, gamma, magnif, kbm, decay;

  double h;

  // Image parameters
  double texp, odttof, det, hfimg, imgbfield;

  // Evap parameters
  double free, p1, t1, tau, beta, p0, offset, t2, tau2, image;
  double finalcpow;
  double trapdepth, hvbar;

  // ODTCALIB parameters
  double odtmaxpow, odtmaxdepth, odtv0radial, odtv0axial;

  // Trap frequencies and geometry factors
  bool w_user;
  double w;
  double wx0, wy0, wz0, a, b, B1, B2, AR;

  // Camera efficiency calibration
  double andoreff_10MHz_14bit_x1_Electron_Mult,
    andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset, eff;

  double mJ_per_photon, photon_per_mJ;

  // Absorption imaging cailbration
  double alphastar;

  // Phase-Contrast analysis
  bool phc;
  double phcsign;

  // Fluorescence analysis
  bool fluor;

  // Analyze pictures from second camera
  bool andor2;

  // Only calculate column density, do not attempt fits
  bool onlyCD;

};

double
b_i (double w, double t)
{
  return w * w / (1 + w * w * t * t);
};

void
init_params (struct params *p)
{
  /************* PARAMETERS OBTAINED FROM REPORT *******/
  p->texp = (double) getINI_num (p->reportfile, "ANDOR", "exp") * 1e-3;	// Exposure in seconds
  p->odttof = (double) getINI_num (p->reportfile, "ODT", "odttof");	// odttof in ms

  p->hfimg = (double) getINI_num (p->reportfile, "ANDOR", "hfimg");	// hfimg in MHz

  double hfimg0 = (double) getINI_num (p->reportfile, "ANDOR", "hfimg0");	// State |1> resonance frequency in MHz
  if (hfimg0 < -1.e10)
    {
      printf
	("Used ANDOR:phcdet to determine phase-contrast detuning, because hfimg0 does not exist\n");
      p->det = (double) getINI_num (p->reportfile, "ANDOR", "phcdet");	// phase contrast detuning in MHz
    }
  else
    {
      p->det = hfimg0 - p->hfimg;
    }




  p->phc = (bool) getINI_num (p->reportfile, "ANDOR", "phc");	// phase contrast 
  p->phcsign = (double) getINI_num (p->reportfile, "ANDOR", "phcsign");	// prefactor in the calculation of the phase contrast signal
  p->fluor = (bool) getINI_num (p->reportfile, "ANDOR", "fluor");	// fluorescence counting

  double Nframes = (double) getINI_num (p->reportfile, "ANDOR", "saveframes");
  if (Nframes < -1.e10)
    {
      printf
	("Used 4 frames for ANDOR.  Could not find ANDOR:saveframes key to determine number of frames\n");
      p->Nframes = 4;
    }
  else
    {
      p->Nframes = (int) floor (Nframes);
    }
  if (p->Nframes == 2)
    {
      p->fluor = true;
    }

  //If pictures come from Andor2, then fluorescence analysis is done
  if (p->andor2)
    {
      p->fluor = true;
      printf ("\n--- ANDOR 2 PICTURES ---\n");
      p->shotnum_fileout = p->shotnum + "_andor2";
    }
  else
    p->shotnum_fileout = p->shotnum;



  //If pictures are blanks then it does not do any fitting
  bool auxblanks = (bool) getINI_num (p->reportfile, "ANDOR", "blanks");	// taking blanks ?
  p->blanks = p->blanks || auxblanks;	// if either the command line or the report says it, these are blanks

//  p->det = -240.;
//  p->phc = 0.;
//  p->phcsign = 1.;
  double magnif_INI = (double) getINI_num (p->reportfile, "ANDOR", "magnif");
  if (magnif_INI == INI_ERROR)
    printf
      ("\n   Did not find ANDOR:magnif in the report file.\n   Using hard-coded value %.2f um/pixel\n\n",
       p->magnif);
  else
    p->magnif = magnif_INI;

  if (p->magnif_override)
    {
      p->magnif = p->magnif_user;
    }

  p->free = (double) getINI_num (p->reportfile, "EVAP", "free");	// 
  p->p1 = (double) getINI_num (p->reportfile, "EVAP", "p1");	// 
  p->t1 = (double) getINI_num (p->reportfile, "EVAP", "t1");	// 
  p->tau = (double) getINI_num (p->reportfile, "EVAP", "tau");	// 
  p->beta = (double) getINI_num (p->reportfile, "EVAP", "beta");	// 
  p->p0 = (double) getINI_num (p->reportfile, "ODT", "odtpow0");	// 
  p->offset = (double) getINI_num (p->reportfile, "EVAP", "offset");	// 
  p->t2 = (double) getINI_num (p->reportfile, "EVAP", "t2");	// 
  p->tau2 = (double) getINI_num (p->reportfile, "EVAP", "tau2");	// 
  p->image = (double) getINI_num (p->reportfile, "EVAP", "image");	//
  p->finalcpow = (double) getINI_num (p->reportfile, "EVAP", "finalcpow");
  p->odtmaxpow = (double) getINI_num (p->reportfile, "ODTCALIB", "maxpow");
  p->odtmaxdepth =
    (double) getINI_num (p->reportfile, "ODTCALIB", "maxdepth");
  p->odtv0radial =
    (double) getINI_num (p->reportfile, "ODTCALIB", "v0radial");
  p->odtv0axial = (double) getINI_num (p->reportfile, "ODTCALIB", "v0axial");


  /************* HARD-CODED PARAMETER VALUES ***********/


  p->eigenface_done = false;

  p->lambda = 670.977e-7;	// cm
  p->hc = 1.98644521e-25;	// J*m 
  p->gamma = 5.9e6;		// Hertz
  p->decay = 27.2e-9;		// Seconds

  p->mJ_per_photon = 1.e3 * (p->hc / (p->lambda * 1.e-2));
  p->photon_per_mJ = 1. / p->mJ_per_photon;


  p->kbm = 1385.;		// um^2/ms^2/uK This is kb/m

  p->probewaist = 0.275;	// cm
  if (p->probewaist_override)
    {
      printf ("override probe waist \n");
      p->probewaist = p->probewaist_user;
    }


  /*** OLD METHOD FOR DETERMINING ANDOR EFFICIENCY - BEFORE SEPTEMBER 2012 ***/
  /***double emp = 0.88;
  p->andoreff_10MHz_14bit_x1_Electron_Mult =
    p->mJ_per_photon * (67.9 / 0.95 / emp);
  // mJ/count = (mJ per photon) * (electrons per A/D count  / QuantumEff / EmpiricalCorrectionFactor) 
  emp = 0.95;
  p->andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset =
    p->mJ_per_photon * (12.9 / 0.95 / emp);***/
  //p->eff = p->andoreff_10MHz_14bit_x1_Electron_Mult;
  //p->eff = p->andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset;


  /*** ANDOR EFFICIENCY CALIBRACTION - SEPTEMBER 2012 ***/
  //  measured probe beam waist =  2750 um 
  //  *Russ Book #4, Page 99. 
  //  *Also in /home/russhart/Probe beam and Camera/ProbeSize 09 06 11.nb)
  //  On 120827 we took pictures at different powers for the probe beam
  //  using a magnification M = 5.  From those pictures (see Pedro Book #5,
  //  Page 232)  we obtained a calibration of the photons/count for the
  //  Andor 
  //  Peak photons incident on camera per pixel per us:
  //  2*P/(pi*(M*w_probe)^2) * (1e-6 sec/usec) * (16 um/pixel)^2 / ( 2.96e-16 mJ/count)  
  //  = 2912.2 * ( P in mW ) 
  //  The peak counts per pixel per us were measured as a function of power
  //  and  were fitted to a line (see excel file in 
  //  /home/russhart/Probe beam and Camera/)
  //  = 111.14 * ( P in mW )  using peak counts
  //  = 142.48 * ( P in mW )  using 10 pixel smoothed data
  //  Putting together the last two above we get
  //  2912.2 / 111.14 = 26.2 photons/count
  //  2912.2 / 142.48 = 20.4 photons/count
  //  From this we obtain the efficiency of the camera:
  //p->eff = pow (0.275 / p->probewaist, 2) * 26.2 * 2.96e-16;  // mJ/count
  p->eff = 26.2 * 2.96e-16;	// mJ/count

  //*** THE MAGNIFICATION SHOULD NOW BE IN ANDOR:magnif ***//
  // 16um/pixel for the Andor, using the 5x obj with the telephoto at f=200
  //p->magnif = 16. / 5.;


  p->h = 48.;			// uK/MHz  this is Planck's constant divided by Boltzmann's constant 
  // Obtain trap depth from report
  p->trapdepth = p->finalcpow / 10. * p->odtmaxdepth;


  // Calculate the geometric mean of trap frequencies
  // Useful to obtain the Fermi temperature
  // EF = h vbar (6N)^(1/3)
  p->hvbar = p->h * pow (p->odtv0radial * p->odtv0radial * p->odtv0axial, 1. / 3.) * sqrt (p->finalcpow / 10.) * 1e-6;	// 1e-6 is to give the trap freqs in MHz
  if (VERBOSE)
    {
      printf
	("..........  DETERMINE GEOMETRIC MEAN OF TRAPPING FREQUENCIES ..........\n");
      printf (" Initial trap depth = %.3f uK\n", p->odtmaxdepth);
      printf ("   Final trap depth = %.3f uK\n", p->trapdepth);
      printf ("   Final trap depth = %.3f %%\n", p->finalcpow / 10. * 100.);
      printf ("  Start radial freq = %.3f Hz\n", p->odtv0radial);
      printf ("   Start axial freq = %.3f Hz\n", p->odtv0axial);
      printf ("  Final radial freq = %.3f Hz\n",
	      p->odtv0radial * sqrt (p->finalcpow / 10.));
      printf ("   Final axial freq = %.3f Hz\n",
	      p->odtv0axial * sqrt (p->finalcpow / 10.));
      printf (" h * ( v_r * v_r * va )^1/3 = %.3f uK\n", p->hvbar);
    }


  //DEBUG: 
  //printf("trap depth = %f\n", p->trapdepth); 
  //printf("hvbar      = %f\n", p->hvbar); 
  // Trap frequencies in kilohertz
  if (!p->w_user)
    p->w = 3.800;

  p->wx0 = p->w * 2 * M_PI;
  p->wy0 = p->w * 2 * M_PI;
  p->wz0 = p->w * 2 * M_PI / 8.;

  p->a = cos (52.5 * M_PI / 180.);
  p->b = sin (52.5 * M_PI / 180.);


  double ws = sqrt (p->finalcpow / 10.);
  double wx = ws * p->wx0;
  double wy = ws * p->wy0;
  double wz = ws * p->wz0;

  //  b(w,t) = w^2/(1+w^2*t^2) 

  double bx = b_i (wx, p->odttof);
  double by = b_i (wy, p->odttof);
  double bz = b_i (wz, p->odttof);

  double a = p->a * p->a;
  double b = p->b * p->b;

  p->B1 =
    -1 * (a * b * (bz - bx) * (bz - bx) / (b * bx + a * bz) - a * bx -
	  b * bz);
  p->B2 = by;


  //DEBUG 
  if (p->show_B)
    {
      //printf("%s: w = %.3f , 1/t^2 = %.5f , B1 = %.5f, B2 = %.5f\n", (p->shotnum).c_str(), p->w, 1/(p->odttof*p->odttof), p->B1, p->B2); 
      printf ("%.3f\t%.5f\t%.5f\t%.5f\n", p->w, 1 / (p->odttof * p->odttof),
	      p->B1, p->B2);
      exit (0);
    }
};


  /********************************************
  PHASE-CONTRAST IMAGING CALCULATION
  ********************************************/

struct phc_calc
{
  double alpha;
  double alpha_pi;
  double phi;
  double phi_pi;
  double atoms;
  double noatoms;
};



class Fermions
{
public:

  Fermions (struct params *params)
  {
    p = params;
  }
   ~Fermions ()
  {
    return;
  }

  void LoadFITS ();
  double ncol_phcimg (double cd, double sig, double det, double i0,
		      double imgbfield, bool twostates, bool SHOW,
		      struct phc_calc *c);
  double signal_phcimg (double ncol, double det, double i0, double imgbfield,
			bool twostates, bool SHOW, struct phc_calc *c);
  void ComputeColumnDensity ();
  void SaveColumnDensity ();
  void FindMoments ();
  void MomentsCrop ();
  void MinimalCrop (double CROP_FACTOR);
  void CropAll (unsigned int roi[4]);
  double GetNPixels ()
  {
    return columndensity->size1 * columndensity->size2;
  }
  void Fit2DGauss (bool mott);
  void FitScatt2DGauss ();
  void FitProbe2DGauss ();
  void Fit2DFermi ();
  void ComputeIntegrated1DDensity ();
  void GetAzimuthalAverageEllipse ();
  void FitAzimuthalFermi ();
  void MakePlots ();

  struct params *p;

  // Quantities obtained in computation of column density 
  double maxI;
  double maxIsmooth;
  double aveI;
  double aveIweighted;
  double maxOD;
  double maxCD;
  double maxPHI;
  double maxSP;
  double maxNSP;
  double maxPEE;
  double maxCA;
  double minCA;
  double maxCN;
  double minCN;
  double Tp0;
  double Tsp;
  double Nsp;
  double maxPHCSIG;		// Maximum phase-contrast signal

  // Arrays for fit results 
  double gaus2dfit[6];
  double gaus2dfit_mott[6];
  double gaus2dfit_err[6];
  double scatt2dfit[6];
  double probe2dfit[5];
  double fermi2dfit[7];
  double fermi_azimuth_fit[5];
  double fermi_azimuth_fit_zero[4];
  double fit1d_gaus_0[4];
  double fit1d_gaus_1[4];
  double fit1d_fermi_0[5];
  double fit1d_fermi_1[5];

  // Quantities obtained as results of fits
  double abs_ci, abs_cj;	// centers of cloud in the uncropped pict

  double number_fit;
  double nfit;
  double nfit_mott;
  double nfit_err;
  double peakd;
  double peakd_mott;
  double peakd_sph;

  double TF;

  double TF_rd, TF_ax, TF_2d, TF_az;	// T/TF for various Fermi fits obtained from BetaMu
  double T_rd, T_ax, T_2d_rd, T_2d_ax, T_az;	// T for various Fermi fits obtained from cloud size


  double n0, BetaMu, cj_Fermi, rj_Fermi, ci_Fermi, ri_Fermi, B_Fermi, Fugacity_Fermi, f_Fermi;	// 2D Fermi parameters 
  double n0_az, BetaMu_az, r_az, B_az, mx_az, Fugacity_az, f_az;	// Azimuthal Fermi parameters
  double n0_az_zeroT, r_az_zeroT, B_az_zeroT, mx_az_zeroT;	// Azimuthal Zero T Fermi parameters

private:
  gsl_matrix * atoms;
  gsl_matrix *noatoms;
  gsl_matrix *atomsref;
  gsl_matrix *noatomsref;

  gsl_matrix *probe;

  gsl_matrix *columndensity;
  gsl_matrix *residuals;

  gsl_matrix *columndensity_scattered_ph;
  gsl_matrix *missing_counts;

  gsl_vector *cutIdata;
  gsl_vector *cutJdata;
  gsl_vector *cutIgauss;
  gsl_vector *cutJgauss;
  gsl_vector *cutIfermi;
  gsl_vector *cutJfermi;

  gsl_vector *sum_density_0;
  gsl_vector *sum_density_1;

  gsl_vector *sum_missing_0;
  gsl_vector *sum_missing_1;

  gsl_vector *sum_density_0_dist;
  gsl_vector *sum_density_1_dist;

  gsl_vector *sum_density_0_fit_gaus;
  gsl_vector *sum_density_1_fit_gaus;

  gsl_vector *sum_density_0_fit_fermi;
  gsl_vector *sum_density_1_fit_fermi;

  unsigned int nbins, usedbins, fitbins;
  double binsize;

  //Arrays for azimuthal averaging
  gsl_vector *azimuthal_all_r;
  gsl_vector *azimuthal_all_dat;

  gsl_vector *azimuthal_r;
  gsl_vector *azimuthal_dat;

  gsl_vector *icut_r;
  gsl_vector *icut_dat;
  gsl_vector *jcut_r;
  gsl_vector *jcut_dat;

  gsl_vector *gaus2d_fit;
  gsl_vector *fermi2d_fit;
  gsl_vector *azimuthal_fit;
  gsl_vector *fermi2d_zero;
  gsl_vector *azimuthal_zero;
  gsl_vector *azimuthal_zero_fit;



  double norm_noat;		//normalization constant for noatoms pict 
  //unsigned int ci, cj, FWHMi, FWHMj, wi1e, wj1e;

};


void
Fermions::LoadFITS ()
{
  /********************************************
  This function loads the FITS data from disk and computes the column density
  matrices.   
 
  It takes care of any necessary cropping as long as it is specified by the user,
  autocropping is  done later.  

  It calls the method ComputeColumnDensity()  to compute the column density.

  ********************************************/

  norm_noat = 1.0;
  abs_ci = 0.;
  abs_cj = 0.;

  if (p->Nframes == 4)
    {
      atoms = ReadFitsImg_gsl_matrix (p->atomsfile);
      noatoms = ReadFitsImg_gsl_matrix (p->noatomsfile);
      atomsref = ReadFitsImg_gsl_matrix (p->atomsreffile);
      noatomsref = ReadFitsImg_gsl_matrix (p->noatomsreffile);
    }
  else if (p->Nframes == 2)
    {
      atoms = ReadFitsImg_gsl_matrix (p->atomsfile);
      noatoms = ReadFitsImg_gsl_matrix (p->noatomsfile);
      atomsref = ReadZeros_gsl_matrix (atoms);
      noatomsref = ReadZeros_gsl_matrix (atoms);
    }
  else
    {
      printf
	(" ERROR:  The program does not know what to do if  Nframes =  %d\n",
	 p->Nframes);
    }



  if (p->crop)
    {
      abs_ci += p->roi[0];
      abs_cj += p->roi[1];

      if (VERBOSE)
	cout << endl << "------------ Cropping Images ------------";
      gsl_matrix *catoms = p->roi_user ? cropImage_ROI (p->roi,
							atoms) :
	cropImage (p->reportfile,
		   atoms);
      gsl_matrix *cnoatoms = p->roi_user ? cropImage_ROI (p->roi,
							  noatoms) :
	cropImage (p->reportfile,
		   noatoms);
      gsl_matrix *catomsref = p->roi_user ? cropImage_ROI (p->roi,
							   atomsref) :
	cropImage (p->reportfile,
		   atomsref);
      gsl_matrix *cnoatomsref = p->roi_user ? cropImage_ROI (p->roi,
							     noatomsref) :
	cropImage (p->reportfile,
		   noatomsref);
      gsl_matrix_free (atoms);
      gsl_matrix_free (noatoms);
      gsl_matrix_free (atomsref);
      gsl_matrix_free (noatomsref);
      atoms = catoms;
      noatoms = cnoatoms;
      atomsref = catomsref;
      noatomsref = cnoatomsref;
    }


  if (VERBOSE)
    cout << endl;

  unsigned int s1 = atoms->size1;
  unsigned int s2 = atoms->size2;

  if (s1 != noatoms->size1 || s2 != noatoms->size2)
    {
      cout << "Atoms and NoAtoms matrix dimensions differ!!" << endl;
      exit (EXIT_FAILURE);
    }

  if (s1 != atomsref->size1 || s2 != atomsref->size2)
    {
      cout << "Atoms and Reference matrix dimensions differ!!" << endl;
      exit (EXIT_FAILURE);
    }

  if (s1 != noatomsref->size1 || s2 != noatomsref->size2)
    {
      cout << "Atoms and noatoms Reference matrix dimensions differ!!" <<
	endl;
      exit (EXIT_FAILURE);
    }

  //Override this values to define a different normalization region
  unsigned int imin, imax = s1;
  unsigned int jmin, jmax = s2;
  imin = imax - imax / 10;
  jmin = jmax - jmax / 10;

  gsl_matrix *norm = gsl_matrix_alloc (imax - imin, jmax - jmin);
  double an = 0.0;		// sum_{i} atoms_{i} * noatoms_{i}
  double nn = 0.0;		// sum_{i} noatoms_{i}^2 

  double a_ = 0.0;		// sum_{i} atoms_{i}
  double n_ = 0.0;		// sum_{i} noatoms_{i} 


  double at, noat, atref, noatref, c0, c1;

  for (unsigned int i = imin; i < imax; i++)
    {
      for (unsigned int j = jmin; j < jmax; j++)
	{

	  gsl_matrix_set (norm, i - imin, j - jmin,
			  gsl_matrix_get (noatoms, i, j));

	  at = gsl_matrix_get (atoms, i, j);
	  noat = gsl_matrix_get (noatoms, i, j);
	  atref = gsl_matrix_get (atomsref, i, j);
	  noatref = gsl_matrix_get (noatomsref, i, j);

	  c0 = noat - noatref;
	  c1 = at - atref;

	  an += c1 * c0;
	  nn += c0 * c0;

	  a_ += c1;
	  n_ += c0;
	}
    }

//  cout << "Norm (ave) = " << a_/n_ << endl; 
//  cout << "Norm (min) = " << an/nn << endl; 

  // Bad idea to use normalization at the moment because noatoms shot still has light 
  // coupled into the fiber for 45 ms when the probe AOM is off
  // norm_noat = an / nn;
  norm_noat = 1.0;

  char base[MAXPATHLEN];
  getcwd (base, MAXPATHLEN);
  string norm_path = makepath (base, p->shotnum_fileout, "_normpatch.TIFF");
  if (p->savetiff)
    toTiffImage (norm, norm_path);
  gsl_matrix_free (norm);

  ComputeColumnDensity ();
  FindMoments ();

  // ROI pos
  double EIGEN_ROI_FACTOR = 3.0;
  double pos0 = max (gaus2dfit[0] - EIGEN_ROI_FACTOR * gaus2dfit[1], 0.);
  double pos1 = max (gaus2dfit[2] - EIGEN_ROI_FACTOR * gaus2dfit[3], 0.);
  double siz0 = min (columndensity->size1 - pos0,
		     gaus2dfit[0] + EIGEN_ROI_FACTOR * gaus2dfit[1] - pos0);
  double siz1 = min (columndensity->size2 - pos1,
		     gaus2dfit[2] + EIGEN_ROI_FACTOR * gaus2dfit[3] - pos1);

  stringstream eigenfacestr;
  eigenfacestr << "eigenface.py ";
  eigenfacestr << p->shotnum << " ";
  eigenfacestr << pos0 << "," << pos1 << "," << siz0 << "," << siz1 << " ";
  eigenfacestr << "100";
  stringstream eigen_background;
  eigen_background << p->shotnum << "_eigenclean.ascii";
  string eigen_filename = eigen_background.str ();

  if (p->eigenface && !p->crop)
    {
      if (VERBOSE)
	{
	  cout << "First few elements of noatoms matrix BEFORE eigenface: " <<
	    endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (noatoms, i, j));
		}
	      cout << "... " << endl;
	    }
	}
      if (VERBOSE)
	{
	  cout << "First few elements of noatomsref matrix BEFORE eigenface: "
	    << endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (noatomsref, i, j));
		}
	      cout << "... " << endl;
	    }
	}


      printf ("%s\n", eigenfacestr.str ().c_str ());
      system (eigenfacestr.str ().c_str ());

      for (unsigned int i = 0; i < s1; i++)
	{
	  for (unsigned int j = 0; j < s2; j++)
	    {
	      gsl_matrix_set (noatomsref, i, j, 0.);
	    }
	}
      noatoms = read_gsl_matrix_ASCII (eigen_filename);
      p->eigenface_done = true;
      ComputeColumnDensity ();
      FindMoments ();

      if (VERBOSE)
	{
	  cout << "First few elements of noatoms matrix AFTER eigenface: " <<
	    endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (noatoms, i, j));
		}
	      cout << "... " << endl;
	    }
	}
      if (VERBOSE)
	{
	  cout << "First few elements of noatomsref matrix AFTER eigenface: "
	    << endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (noatomsref, i, j));
		}
	      cout << "... " << endl;
	    }
	}
    }

  return;
}



  /********************************************
  PHASE-CONTRAST IMAGING CALCULATION
  ********************************************/


struct phc_params
{
  double sig;
  double det;
  double i0;
  double imgbfield;
  bool twostates;
  bool SHOW;
  double lambda;
  double gamma;
  double magnif;
  struct phc_calc *c;
};

//Formula to calculate splitting between states 1 and 2
double
split12 (double bfield)
{
  return 66.718 - 0.009 * bfield + 0.601 * sqrt (bfield);
}


double
sig_phcimg (double ncol, struct phc_params *phc)
{

  //This is the maximal cross section ( units of pixel )
  double sigma0 =
    6 * M_PI * pow (phc->lambda / 2. / M_PI, 2) / pow (phc->magnif * 1e-4, 2);

  //This is the cross section for light polarized perp
  //to the magnetic field drving sigma minus transitions
  //( units of pixel )
  double sig_minus =
    3 * M_PI * pow (phc->lambda / 2. / M_PI, 2) / pow (phc->magnif * 1e-4, 2);

  //This is the cross section for light polarized along
  //to the magnetic field drving deltaM=0 transitions
  //( units of pixel )
  double sig_pi =
    4 * M_PI * pow (phc->lambda / 2. / M_PI, 2) / pow (phc->magnif * 1e-4, 2);

  //This is the cross section for light polarized perp
  //to the magnetic field drving sigma plus transitions
  //( units of pixel )
  double sig_plus =
    1 * M_PI * pow (phc->lambda / 2. / M_PI, 2) / pow (phc->magnif * 1e-4, 2);


  //The I/Isat parameter is obtained for each transition:
  //i0 is I/Isat with Isat = 5.1, i.e. maximal cross section.
  double i0_minus = phc->i0 * sig_minus / sigma0;
  double i0_pi = phc->i0 * sig_pi / sigma0;
  double i0_plus = phc->i0 * sig_plus / sigma0;

  //The detunings are calculated:
  // det is detuning from state |1> for the sigma minus transition
  double det = phc->det;

  // splitting between states |1> and |2>:
  double delta12 = split12 (phc->imgbfield) / (phc->gamma / 1.e6);

  // detuning from state2
  double det2 = det + delta12;

  // detuning for pi transition, state |1>
  double det_pi = det - 1.87 * phc->imgbfield / (phc->gamma / 1.e6);
  // detuning for pi transition, state |2>
  double det2_pi = det_pi + delta12;
  // detuning for plus transition, state |1>
  double det_plus = det - 2 * 1.87 * phc->imgbfield / (phc->gamma / 1.e6);
  // detuning for plus transition, state |2>
  double det2_plus = det_plus + delta12;
  //Here, the parameters for the phase contrast imaging setup are set
  //First the polarization of the probe light
  double a = 1. / sqrt (2.);
  double b = 1. / sqrt (2.);
  double g = M_PI / 2.;
  //Then the angle of the polarizer wrt the magnetic field
  double th = -M_PI / 4.;
  //One atom cannot scatter a photon twice, so it is necessary to 
  //include in the absorption (alpha) and phase-shift (phi) terms
  //a term corresponding to the probability of the atom to undergo
  //each transition
  //Then the absorption and phase shifts are calculated 
  //The unsubscripted one is for the sigma minus transition
  //The pi subscript one is for the pi transition
  //double alpha, alpha_pi, phi, phi_pi;
  if (phc->twostates)
    {
      //Calculate the cross section to scatter on each transition from each state:
      double minus1 = sig_minus / (1 + 4 * det * det + 2 * i0_minus);
      double minus2 = sig_minus / (1 + 4 * det2 * det2 + 2 * i0_minus);
      double pi1 = sig_pi / (1 + 4 * det_pi * det_pi + 2 * i0_pi);
      double pi2 = sig_pi / (1 + 4 * det2_pi * det2_pi + 2 * i0_pi);
      double plus1 = sig_plus / (1 + 4 * det_plus * det_plus + 2 * i0_plus);
      double plus2 = sig_plus / (1 + 4 * det2_plus * det2_plus + 2 * i0_plus);
      //One atom cannot scatter a photon twice, so it is necessary to 
      //include in the absorption (alpha) and phase-shift (phi) terms
      //a term corresponding to the probability of the atom to undergo
      //each transition
      //Calculate the probabilities to scatter on each transition:
      /*
         // If atom is in state |1>:
         double p_minus1 = minus1 / (minus1 + pi1);
         double p_pi1 = pi1 / (minus1 + pi1);
         // If atom is in state |2>:
         double p_minus2 = minus2 / (minus2 + pi2);
         double p_pi2 = pi2 / (minus2 + pi2);
       */
      // If atom is in state |1>:
      /*
         double p_minus1 = sig_minus / (sig_minus + sig_pi);
         double p_pi1 = sig_pi / (sig_minus + sig_pi);
         // If atom is in state |2>:
         double p_minus2 = sig_minus / (sig_minus + sig_pi);
         double p_pi2 = sig_pi / (sig_minus + sig_pi);
       */
      // Calculate the absorption due to atoms in state 1 
      double alpha1 = (ncol / 2.) * minus1;
      double alpha1_pi = (ncol / 2.) * pi1;
      double alpha1_plus = (ncol / 2.) * plus1;
      // and in state 2
      double alpha2 = (ncol / 2.) * minus2;
      double alpha2_pi = (ncol / 2.) * pi2;
      double alpha2_plus = (ncol / 2.) * plus2;
      /* 
         double alpha1 = 
         (ncol / 2.) * sig_minus / (1 + 4 * det * det + 2 * i0_minus);
         double alpha2 =
         (ncol / 2.) * sig_minus / (1 + 4 * det2 * det2 + 2 * i0_minus);

         double alpha1_pi =
         (ncol / 2.) * sig_pi / (1 + 4 * det_pi * det_pi + 2 * i0_pi);
         double alpha2_pi =
         (ncol / 2.) * sig_pi / (1 + 4 * det2_pi * det2_pi + 2 * i0_pi);
       */
      phc->c->alpha = alpha1 + alpha2 + alpha1_plus + alpha2_plus;
      phc->c->alpha_pi = alpha1_pi + alpha2_pi;
      phc->c->phi =
	-alpha1 * det - alpha2 * det2 - alpha1_plus * det_plus -
	alpha2_plus * det2_plus;
      phc->c->phi_pi = -alpha1_pi * det_pi - alpha2_pi * det2_pi;
    }

  else
    {
      double alpha_plus =
	ncol * sig_plus / (1 + 4 * det_plus * det_plus + 2 * i0_plus);
      double alpha_minus =
	ncol * sig_minus / (1 + 4 * det * det + 2 * i0_minus);
      phc->c->alpha = alpha_plus + alpha_minus;
      phc->c->alpha_pi =
	ncol * sig_pi / (1 + 4 * det_pi * det_pi + 2 * i0_pi);
      phc->c->phi = -alpha_minus * det - alpha_plus * det_plus;
      phc->c->phi_pi = -phc->c->alpha_pi * det_pi;
    }

  phc->c->atoms = b * b * exp (-phc->c->alpha_pi) * cos (th) * cos (th)
    + a * a * exp (-phc->c->alpha) * sin (th) * sin (th)
    + a * b * exp (-phc->c->alpha / 2. - phc->c->alpha_pi / 2.) * cos (g -
								       phc->
								       c->
								       phi +
								       phc->
								       c->
								       phi_pi)
    * sin (2. * th);
  phc->c->noatoms =
    b * b * cos (th) * cos (th) + a * a * sin (th) * sin (th) +
    a * b * cos (g) * sin (2 * th);
  return -1. * (phc->c->atoms / phc->c->noatoms - 1);
};

double
sig_phcimg_err (double ncol, void *params)
{
  struct phc_params *phc = (struct phc_params *) params;
  return sig_phcimg (ncol, phc) - phc->sig;
};

double
Fermions::signal_phcimg (double ncol, double det, double i0,
			 double imgbfield, bool twostates, bool SHOW,
			 struct phc_calc *c)
{

  struct phc_params phc = {
    0.0, det, i0, p->imgbfield, true, false, p->lambda, p->gamma,
    p->magnif, c
  };
  return sig_phcimg (ncol, &phc);
}

double
Fermions::ncol_phcimg (double cd, double sig, double det, double i0,
		       double imgbfield, bool twostates, bool SHOW,
		       struct phc_calc *c)
{


  //This is the maximal cross section ( units of pixel )
  //double sigma0 =
  //  6 * M_PI * pow (p->lambda / 2. / M_PI, 2) / pow (p->magnif * 1e-4, 2);

  //First use column density from the linearized algortihm to 
  //make a first estimate of the column density 
  /*
     double ncol = 1.5 * pow (fabs (det / 25),
     1. / 3.) * sig / (0.5 - det) * (1 + 4 * det * det +
     2 * i0) / sigma0; */
  double ncol = 1. * cd * (det < -25 ? pow (fabs (det / 25), 1. / 3.) : 1.);
  if (SHOW)
    {
      cout << endl <<
	"-------- CALCULATING PHASE-CONTRAST COLUMN DENSITY -----" << endl;
      printf (" First naive estimate = %.3f\n", ncol);
    }

  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  double x_lo, x_hi;
  if (fabs (ncol) < 100.)
    {
      x_lo = -1000.;
      x_hi = 1000.;
    }
  else
    {
      double span = 1.0;
      x_lo = ncol - span * fabs (ncol);
      x_hi = ncol + span * fabs (ncol);
    }

  gsl_function F;
  F.function = &sig_phcimg_err;
  struct phc_params params = {
    sig, det, i0, imgbfield, twostates, false, p->lambda, p->gamma,
    p->magnif, c
  };
  F.params = &params;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  if (SHOW)
    {
      printf ("using %s method\n", gsl_root_fsolver_name (s));
      printf ("%5s [%9s, %9s] %9s %9s\n",
	      "iter", "lower", "upper", "root", "err(est)");
    }
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
      if (status == GSL_SUCCESS && SHOW)
	printf ("Converged:\n");
      if (SHOW)
	{
	  printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
		  iter, x_lo, x_hi, r, x_hi - x_lo);
	}
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  double ratio = r / ncol;
  if (fabs (ncol) > 100. && (ratio > 1.7 || ratio < 0.3))
    {
      printf (" Result / FirstEstimate = %.3f\n", ratio);
    }

  return r;
}

void
Fermions::ComputeColumnDensity ()
{
  /********************************************
  Four matrices are calculated here  

  columndensity = the number of atoms in each pixel

  columndensity_scattered_ph = the number of atoms in each pixel, 
                               calculated from the number of photons 
                               that were scattered

  missing_counts = counts in atoms picture minus counts in no-atoms picture

  probe = the number of photons that were incident on each pixel 

  ********************************************/

  unsigned int s1 = atoms->size1;
  unsigned int s2 = atoms->size2;
  columndensity = gsl_matrix_alloc (s1, s2);
  columndensity_scattered_ph = gsl_matrix_alloc (s1, s2);
  missing_counts = gsl_matrix_alloc (s1, s2);
  probe = gsl_matrix_alloc (s1, s2);
  /************ ABSORPTION IMAGING VARIABLES ************/
  double at, noat, atref, noatref, c0, c1, cd, OD, i0, i0eff, p0, term1,
    term2;
  double cd_sp, cd_highintabs;
  double maxT1, maxT2;		//maxima for term1 and term2 in the column density
  double maxCDterm1 = 0., maxCDterm2 = 0.;
  maxT1 = -1e6;
  maxT2 = -1e6;
  double sigma0 = 3 * M_PI * pow (p->lambda / 2 / M_PI, 2);	//This corresponds to half the maximal cross-section, because of our polarization.  
  double isat = 10.2;		// This corresponds to twice the minimal Isat because of our polarization
  //should be equal to  gamma h v / sigma0 .... test this.  
  double det = 0.;		//for now all  absorption imaging is on resonance 
  maxI = 0.;
  aveI = 0.;
  aveIweighted = 0.;
  double sumcd = 0.;
  maxCD = 0.;			//maximum column density
  maxPHI = 0.;			//maximum phase shift in phase constrast
  maxNSP = 0.;
  maxSP = 0.;
  Tp0 = 0.;			//total photons in noatoms frame 
  double sp = 0.;		//scattered photons in pixel 
  Tsp = 0.;			//total number of scattered photons 
  double nsp = 0.;		//number from scattered photons in pixel
  Nsp = 0.;			//total number from scattered photons
  double pee = 0.;
  maxPEE = 0.;
  double avg_pee = 0.;		//average rho_ee
  maxOD = 0.;
  double OD_err = 5;		//any higher OD will be cutoff and trigger a warning
  bool OD_errmsg = false;
  bool high_phase_shift_flag = false;
  //bool sanity_check_flag = false;
  //sanity_check_flag = false;
  minCA = 1.e6;			// minimum counts in atoms frame
  minCN = 1.e6;			// minimum counts in noatoms frame
  maxCA = 0.;			// maximum counts in atoms frame
  maxCN = 0.;			// maximum counts in noatoms frame
  //Structure to store intermediate results of phc calc
  struct phc_calc calc = {
    0., 0., 0., 0., 0., 0.
  };
  struct phc_calc calcMAX = {
    0., 0., 0., 0., 0., 0.
  };
  //For max column density, counts in atoms and noatoms frame
  double c0CDMAX = 0.0;
  double c1CDMAX = 0.0;
  //Matrix for OD
  gsl_matrix *od_matrix;
  od_matrix = gsl_matrix_alloc (s1, s2);

  //Loss of intensity at the polarizer
  double loss = 1.0;
  if (p->phc)
    loss = 0.11 / 0.26;		//Measured on 10/17/2012 Russ #8 Page98

  if (p->cdsp && p->highintabs)
    {
      printf
	("Cannot use --cdsp and --highintabs options at the same time\n");
      printf ("Program will exit\n");
      exit (EXIT_FAILURE);
    }

  printf ("---  CAMERA PARAMETERS: ---\n");
  //printf ("  efficiency    = %.3e mJ/count\n", p->eff);
  //printf ("  photons/mJ    = %.3e\n", p->photon_per_mJ);
  printf ("  frames saved  = %d\n", p->Nframes);
  printf ("  photons/count = %.3e\n", p->eff * p->photon_per_mJ);
  printf ("  magnif        = %.2f um/pixel\n", p->magnif);

  double probeintensity;

  /************ PHASE-CONTRAST IMAGING VARIABLES ************/
  double signal, phase_shift, sig1, sig2, state2det, state1det, delta12;
  gsl_matrix *phc_signal;
  if (!p->phc && !p->fluor)
    {
      printf ("---  ABSORPTION IMAGING ---\n");
      if (p->cdsp)
	{
	  printf
	    ("... User selected Col Dens. From Scattered Photons calculation.\n");
	}
      if (p->highintabs)
	{
	  printf
	    ("... User selected High-Intensity Absorption imaging calculation.\n");
	  printf ("... Below are the probe parameters that will be used. \n");
	  probeintensity =
	    2. * p->probepower / (M_PI * pow (p->probewaist, 2));
	  printf ("  probepower    = %.3f mW\n", p->probepower);
	  printf ("  probewaist    = %.3f cm\n", p->probewaist);
	  printf ("  probeintensity= %.3f isat\n", probeintensity / isat);
	}
    }
  else if (p->fluor)
    {
      printf ("--- FLUORESCENCE IMAGING ---\n");
    }
  else if (p->phc)
    {
      phc_signal = gsl_matrix_alloc (s1, s2);	// A matrix of the phase-contrast Signal at each pixel.
      sigma0 = 3 * M_PI * pow (p->lambda / 2. / M_PI, 2);
      isat = 10.2;
      det = p->det * 1.e6 / p->gamma;	//detuning in units of gamma
      state1det = det;
      p->imgbfield = (-1. * (p->hfimg + p->det) - 100.0 - 163.7) / -1.414;
      delta12 = split12 (p->imgbfield) * 1.e6 / p->gamma;	// Splitting between states |1> and |2> 
      state2det = state1det + delta12;

      printf ("---  PHASE-CONTRAST PARAMETERS: ---\n");
      if (p->twostates)
	printf ("  two states = true\n");
      else
	printf ("  two states = false\n");
      printf ("  img bfield = %.2f Gauss\n", p->imgbfield);
      printf ("  hfimg      = %.2f MHz\n", p->hfimg);
      printf ("  Phc det    = %+-3.2f Gamma  = %+-3.2f MHz\n",
	      det, det * p->gamma / 1.e6);
      printf ("  |1> to |2> = %+-3.2f Gamma  = %+-3.2f MHz\n",
	      delta12, delta12 * p->gamma / 1.e6);
    }


  /********************************************
  *
  *  Iterate over the image matrices to calculate the column
  *  density in each pixel
  *
  *********************************************/
  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{

	  at = gsl_matrix_get (atoms, i, j);
	  noat = gsl_matrix_get (noatoms, i, j);
	  atref = gsl_matrix_get (atomsref, i, j);
	  noatref = gsl_matrix_get (noatomsref, i, j);
	  if (p->eigenface && noatref != 0. && p->eigenface_done)
	    {
	      printf ("eigenface error: noatref != 0, noatref = %.3f\n",
		      noatref);
	    }

	  if (p->Nframes == 2)
	    {
	      // Ensure theres is no probe normalization
	      // Ensure reference frames are considered zero
	      norm_noat = 1.;
	      noatref = 0.;
	      atref = 0.;
	    }

	  c0 = norm_noat * (noat - noatref);	// counts in no atoms frame
	  c1 = at - atref;	// counts in atoms frame

	  if (c0 < minCN)
	    minCN = c0;
	  if (c0 > maxCN)
	    maxCN = c0;
	  if (c1 < minCA)
	    minCA = c1;
	  if (c1 > maxCA)
	    maxCA = c1;
	  // eff is in mJ/count
	  p0 = c0 * p->eff * p->photon_per_mJ;	// number of photons in pixel in noatoms frame
	  Tp0 += p0;		// total number of photons in noatoms frame 



	  i0 = (1 / loss) * (c0 * p->eff / (pow (p->magnif * 1e-4, 2)) / p->texp / (isat));	// intensity in noatoms frame in units of isat
	  i0eff = i0 / p->alphastar;	// effective intensity in noatoms frame

	  aveI += i0 / s1 / s2;
	  if (i0 > maxI)
	    {
	      maxI = i0;
	      double Ismooth = 0.;
	      int nsmooth = 0;
	      int dsmooth = 10;
	      for (int ii = -dsmooth; ii <= dsmooth; ii++)
		{
		  for (int jj = -dsmooth; jj <= dsmooth; jj++)
		    {
		      if (i + ii >= 0 && i + ii < s1 && j + jj >= 0
			  && j + jj < s2)
			{
			  double c0temp = gsl_matrix_get (noatoms, i + ii,
							  j + jj) -
			    gsl_matrix_get (noatomsref, i + ii, j + jj);
			  Ismooth +=
			    (c0temp * p->eff / (pow (p->magnif * 1e-4, 2)) /
			     p->texp / (isat));
			  nsmooth++;
			}
		    }
		}
	      maxIsmooth = Ismooth / nsmooth;
	    }


	  OD = log (fabs (c0 / c1));
	  gsl_matrix_set (od_matrix, i, j, OD);
	  if (OD > OD_err and ! p->fluor)
	    {
	      OD = OD_err;
	      if (!OD_errmsg)
		printf
		  ("\n******   WARNING: Optical density too high !!! ******\n");
	      OD_errmsg = true;
	    }

	  else if (OD > maxOD)
	    maxOD = OD;


	  /*********** ABSORPTION IMAGING ***************/
	  if (!p->phc && !p->fluor)
	    {

	      term1 =
		p->alphastar * (1. +
				4. * det * det) * OD / sigma0 *
		pow (p->magnif * 1e-4, 2);
	      pee = i0eff / (1 + 2 * i0eff);
	      avg_pee += pee / (s1 * s2);
	      if (pee > maxPEE)
		maxPEE = pee;
	      //sp = scattered photons in pixel
	      sp = (c0 - c1) * p->eff * p->photon_per_mJ;
	      //sp = (c0 - c1) * 24.2;
	      if (sp > maxSP)
		maxSP = sp;
	      //nsp = number of atoms in pixel, calculated from number of scattered photons
	      nsp = sp / (pee * p->texp / p->decay);
	      if (nsp > maxNSP)
		maxNSP = nsp;
	      Tsp += sp;	// total number of scattered photons
	      Nsp += nsp;	// total number of atoms from scattered photons 
	      term2 =
		2 * i0 * (1. - c1 / c0) / sigma0 * pow (p->magnif * 1e-4, 2);
	      if (term1 > maxT1)
		maxT1 = term1;
	      if (term2 > maxT2)
		maxT2 = term2;
	      cd = (term1 + term2);

	      //Calculate column density from scattered photons
	      cd_sp =
		(1 -
		 c1 / c0) * ((1 +
			      2 * i0) / i0) * (p->decay /
					       p->mJ_per_photon) *
		pow (p->magnif * 1e-4, 2) * i0 * isat;

	      //Calculate column density from scattered photons in the high intensity approximation
	      //with probe power and beam waist provided by user
	      cd_highintabs =
		(1 -
		 c1 / c0) * (2. * p->decay / p->mJ_per_photon) *
		pow (p->magnif * 1e-4, 2) * (probeintensity);
	      if (p->cdsp)
		{
		  cd = cd_sp;
		}
	      if (p->highintabs)
		{
		  cd = cd_highintabs;
		}
	    }

	  /*********** FLUORESCENCE IMAGING ***************/
	  else if (p->fluor)
	    {

	      cd = c1 - c0;
	      nsp = 0.;		//this is not relevant for phase contrast so just make it zero.

	    }

	  /*********** PHASE - CONTRAST IMAGING ***************/
	  else
	    {
	      signal = p->phcsign * (c1 / c0 - 1);
	      gsl_matrix_set (phc_signal, i, j, signal);
	      phase_shift = asin (signal);
	      //phase_shift = -1.0 * det * signal / (0.5 - det);
	      //cout << "signal = " << signal << " phi = " << phase_shift << endl; 
	      if (fabs (phase_shift) > maxPHI)
		maxPHI = phase_shift;
	      //if (signal / (0.5 - det) > 0.5 || signal > 0.6)
	      if (phase_shift > M_PI / 5.)
		{
		  if (!high_phase_shift_flag)
		    printf
		      ("\n******   WARNING: Maximum phase shift exceeds pi/5 !!! ******\n");
		  high_phase_shift_flag = true;
		}



	      /* THIS IS THE OLD PHASE-CONTRAST CALCULATION */
	      sig1 =
		(0.5 - state1det) * sigma0 / pow (p->magnif * 1e-4,
						  2) / (1 +
							4 * state1det *
							state1det + 2 * i0);
	      sig2 =
		(0.5 - state2det) * sigma0 / pow (p->magnif * 1e-4,
						  2) / (1 +
							4 * state2det *
							state2det + 2 * i0);
	      if (p->twostates)
		cd = 2 * signal / (sig1 + sig2);	//columdensity
	      else
		cd = 2 * signal / (sig1);	//columdensity -> IF USING ONLY STATE 1 

	      /* THIS IS THE NEW PHASE-CONTRAST CALCULATION */
	      bool SHOW = false;
	      cd =
		ncol_phcimg (cd, signal, det, i0, p->imgbfield, p->twostates,
			     SHOW, &calc);
	      if (i == 302 && j == 200 && false)
		{
		  printf
		    ("\nDEBUG PHASE-CONTRAST COLUMN DENSITY CALCULTION\n");
		  printf (" i = %u, j = %u\n", i, j);
		  printf ("    signal = %.6f\n", signal);
		  printf ("        cd = %.3f\n", cd);
		  printf ("signal(cd) = %.6f\n",
			  signal_phcimg (cd, det, i0, p->imgbfield, true,
					 true, &calc));
		}


	      nsp = 0.;		//this is not relevant for phase contrast so just make it zero.
	      // sanity check
	      /*if ( !sanity_check_flag){
	         printf("det    = %.3e * Gamma\n", det );
	         printf("det^2  = %.3e * Gamma^2\n", det*det );
	         printf("i0     = %.3e * Isat\n", i0 );
	         printf("sigma0 = %.3e\n", sigma0 );
	         printf("sig1   = %.3e\n", sig1);
	         printf("sig2   = %.3e\n", sig2);
	         printf("cd     = %.3e\n", cd);
	         printf("signal = %.3e\n",signal);
	         } */
	    }
	  if (cd > maxCD)
	    {
	      maxCD = cd;
	      c0CDMAX = c0;
	      c1CDMAX = c1;
	      if (p->phc)
		calcMAX = calc;
	      if (!p->phc && !p->fluor)
		{
		  maxCDterm1 = term1;
		  maxCDterm2 = term2;
		}
	    }

	  if (cd > 0.)
	    {
	      aveIweighted += i0 * cd;
	      sumcd += cd;
	    }

	  gsl_matrix_set (columndensity_scattered_ph, i, j, nsp);
	  gsl_matrix_set (columndensity, i, j, cd);
	  gsl_matrix_set (missing_counts, i, j, c1 - c0);
	  gsl_matrix_set (probe, i, j, c0);
	}
    }

  aveIweighted = aveIweighted / sumcd;

  if (!p->fluor)
    {
      unsigned int smooth_bins = 4;
      gsl_matrix *od_matrix_smoothed = smooth (od_matrix, smooth_bins);
      unsigned int odmax_i, odmax_j;
      double od_max_pos;
      findpeak (od_matrix_smoothed, &odmax_i, &odmax_j, &od_max_pos, true);
      maxOD = od_max_pos;
    }
  if (p->phc && !p->fluor)
    {
      unsigned int smooth_bins = 4;
      gsl_matrix *phc_signal_smoothed = smooth (phc_signal, smooth_bins);
      unsigned int phcsig_i, phcsig_j;
      double phcsig_max_pos, phcsig_max_neg;
      findpeak (phc_signal_smoothed, &phcsig_i, &phcsig_j,
		&phcsig_max_pos, true);
      findpeak (phc_signal_smoothed, &phcsig_i, &phcsig_j,
		&phcsig_max_neg, false);
      maxPHCSIG =
	fabs (phcsig_max_pos) >
	fabs (phcsig_max_neg) ? phcsig_max_pos : phcsig_max_neg;
    }



  if (VERBOSE || p->imgverbose)
    {
      cout << endl << "------------ Column Density Stats ------------" <<
	endl;
      if (!p->phc && !p->fluor)
	{
	  cout << endl << "----> Method used : absorption" << endl;
	  printf ("\tSample contains: %s\n",
		  p->twostates ? "spin mixture of |1> and |2>" :
		  "only atoms in state |1>");
	  cout << "\ts1=" << s1 << ", s2=" << s2 << endl;
	  cout << "\talphastar = " << p->alphastar << endl;
	  cout << "\tmax term1 (OD) = " << maxT1 << endl;
	  cout << "\tmax term2 (sp) = " << maxT2 << endl;
	  cout << "\tmax OD = " << maxOD << endl;
	  cout << "\tmax probe intensity = " << maxI << " Isat " << endl;
	  printf ("\tmax probe intensity = %.3f mW/cm^2\n", maxI * isat);
	  printf ("\tmax CD = %.5f ( = %.3f + %.3f)\n", maxCD, maxCDterm1,
		  maxCDterm2);
	  printf ("\tmax CD, cNOATOMS = %.3f\n", c0CDMAX);
	  printf ("\tmax CD, cATOMS = %.3f\n", c1CDMAX);
	  cout << "\tmax scattered photons = " << maxSP << endl;
	  cout << "\tmax number from scattered photons = " << maxNSP << endl;
	  cout << "\tmax rho_{ee} = " << maxPEE << endl;
	  cout << "\tavg rho_{ee} = " << avg_pee << endl;
	  printf ("\tcounts in atoms picture:     ( %f to %f )\n",
		  minCA, maxCA);
	  printf ("\tcounts in no atoms picture:  ( %f to %f )\n",
		  minCN, maxCN);
	  cout << "\ttotal scattered photons = " << Tsp << endl;
	  cout << "\ttotal number from scattered photons = " << Nsp << endl;
	  cout << "\ttotal number of photons in probe pulse = " << Tp0
	    << endl;
	  if (VERBOSE)
	    {
	      cout << "First few elements of column density matrix: " << endl;
	      for (int i = 0; i < 5; i++)
		{
		  for (int j = 0; j < 5; j++)
		    {
		      printf ("%.5f\t ",
			      gsl_matrix_get (columndensity, i, j));
		    }
		  cout << "... " << endl;
		}
	    }

	}
      else if (p->fluor)
	{
	  cout << endl << "----> Method used : flourescence counting photons"
	    << endl;
	  printf ("\tmax photons in pixel = %.5f \n", maxCD);
	  printf ("\tcounts in atoms picture:     ( %f to %f )\n",
		  minCA, maxCA);
	  printf ("\tcounts in no atoms picture:  ( %f to %f )\n",
		  minCN, maxCN);
	}
      else
	{
	  cout << endl << "----> Method used : phase contrast" << endl;
	  printf ("\tSample contains: %s\n",
		  p->twostates ? "spin mixture of |1> and |2>" :
		  "only atoms in state |1>");
	  cout << "\ts1=" << s1 << ", s2=" << s2 << endl;
	  cout << "\tmax probe intensity = " << maxI << " Isat " << endl;
	  printf ("\tmax probe intensity = %.3f mW/cm^2\n", maxI * isat);
	  cout << "\tmax PHI = " << maxPHI << endl;
	  cout << "\tmax CD = " << maxCD << endl;
	  printf ("\tmax CD, cATOMS = %.3f\n", c1CDMAX);
	  printf ("\tmax CD, cNOATOMS = %.3f\n", c0CDMAX);
	  printf ("\tmax CD, alpha    = %.6f\n", calcMAX.alpha);
	  printf ("\tmax CD, alpha_pi = %.6f\n", calcMAX.alpha_pi);
	  printf ("\tmax CD, phi      = %.6f\n", calcMAX.phi);
	  printf ("\tmax CD, phi_pi   = %.6f\n", calcMAX.phi_pi);
	  printf ("\tmax CD, atoms    = %.6f\n", calcMAX.atoms);
	  printf ("\tmax CD, noatoms  = %.6f\n", calcMAX.noatoms);
	  printf ("\tcounts in atoms picture:     ( %f to %f )\n",
		  minCA, maxCA);
	  printf ("\tcounts in no atoms picture:  ( %f to %f )\n",
		  minCN, maxCN);
	  cout << "\ttotal number of photons in probe pulse = " << Tp0
	    << endl;
	  if (VERBOSE)
	    {
	      cout << "First few elements of column density matrix: " << endl;
	      for (int i = 0; i < 5; i++)
		{
		  for (int j = 0; j < 5; j++)
		    {
		      printf ("%.5f\t ",
			      gsl_matrix_get (columndensity, i, j));
		    }
		  cout << "... " << endl;
		}
	    }

	}

    }
  return;
}

void
Fermions::SaveColumnDensity ()
{
  /********************************************
  Four matrices are saved here 

  ... The same four that are created inside Fermions::ComputeColumnDensity() 

  ********************************************/
  char base[MAXPATHLEN];
  getcwd (base, MAXPATHLEN);
  string column_path = makepath (base, p->shotnum_fileout, "_column.TIFF");
  string scatt_path =
    makepath (base, p->shotnum_fileout, "_column_scatt.TIFF");
  string missing_path =
    makepath (base, p->shotnum_fileout, "_missing_counts.TIFF");
  string probe_path = makepath (base, p->shotnum_fileout, "_probe.TIFF");
  string column_ascii_path =
    makepath (base, p->shotnum_fileout, "_column.ascii");
  string scatt_ascii_path =
    makepath (base, p->shotnum_fileout, "_column_scatt.ascii");
  string missing_ascii_path =
    makepath (base, p->shotnum_fileout, "_missing_counts.ascii");
  string probe_ascii_path =
    makepath (base, p->shotnum_fileout, "_probe.ascii");
  if (p->savetiff)
    {
      toTiffImage (columndensity, column_path);
      toTiffImage (columndensity_scattered_ph, scatt_path);
      toTiffImage (probe, probe_path);
      toTiffImage (missing_counts, missing_path);
    }

  save_gsl_matrix_ASCII (columndensity, column_ascii_path);
  if (p->saveascii)
    {
      save_gsl_matrix_ASCII (columndensity_scattered_ph, scatt_ascii_path);
      save_gsl_matrix_ASCII (probe, probe_ascii_path);
      save_gsl_matrix_ASCII (missing_counts, missing_ascii_path);
    }
  return;
}



void
Fermions::CropAll (unsigned int roi[4])
{
  gsl_matrix *cropped_columndensity = cropImage_ROI (roi, columndensity);
  gsl_matrix *cropped_columndensity_scattered_ph =
    cropImage_ROI (roi, columndensity_scattered_ph);
  gsl_matrix *cropped_probe = cropImage_ROI (roi, probe);
  gsl_matrix *cropped_missing_counts = cropImage_ROI (roi, missing_counts);
  gsl_matrix_free (columndensity);
  gsl_matrix_free (columndensity_scattered_ph);
  gsl_matrix_free (probe);
  gsl_matrix_free (missing_counts);
  columndensity = cropped_columndensity;
  columndensity_scattered_ph = cropped_columndensity_scattered_ph;
  probe = cropped_probe;
  missing_counts = cropped_missing_counts;
  return;
}

void
Fermions::FindMoments ()
{
  Gaus2DGuess (columndensity, gaus2dfit, p->shotnum_fileout, false);
  return;
}


void
Fermions::MinimalCrop (double CROP_FACTOR)
{
  /********************************************
  Three matrices are cropped here 
 
  ... The same four that are created inside Fermions::ComputeColumnDensity() 

  ********************************************/
  if (VERBOSE)
    cout << endl <<
      "------------ CROPPING COLUMN DENSITY TO MINIMAL AREA------------"
      << endl;
  // Results from 2DGaus Fit are used to crop the column density
  // center is (ci_,cj_)  
  // sizes are (wi_1e, wj_1e)
  // ROI is defined as (ax0pos, ax1pos, ax0size, ax1size)
  unsigned int roi[4];
  // ROI pos
  double pos0 = max (gaus2dfit[0] - CROP_FACTOR * gaus2dfit[1], 0.);
  double pos1 = max (gaus2dfit[2] - CROP_FACTOR * gaus2dfit[3], 0.);
  double siz0 = min (columndensity->size1 - pos0,
		     gaus2dfit[0] + CROP_FACTOR * gaus2dfit[1] - pos0);
  double siz1 = min (columndensity->size2 - pos1,
		     gaus2dfit[2] + CROP_FACTOR * gaus2dfit[3] - pos1);
  roi[0] = (unsigned int) floor (pos0);
  roi[1] = (unsigned int) floor (pos1);
  roi[2] = (unsigned int) floor (siz0);
  roi[3] = (unsigned int) floor (siz1);
  if (VERBOSE)
    {
      printf
	("    Determined ROI for minimal crop [%d,%d,%d,%d]\n",
	 roi[0], roi[1], roi[2], roi[3]);
    }

  abs_ci += roi[0];
  abs_cj += roi[1];
  gaus2dfit[0] = gaus2dfit[0] - double (roi[0]);
  gaus2dfit[2] = gaus2dfit[2] - double (roi[1]);
  CropAll (roi);
  if (VERBOSE)
    {
      printf ("\n    New matrix dimensions = %d, %d\n\n",
	      (unsigned int) columndensity->size1,
	      (unsigned int) columndensity->size2);
    }
  return;
}

void
Fermions::Fit2DGauss (bool mott = 0)
{

  if (VERBOSE)
    cout << endl <<
      "------------ FIT COLUMN DENSITY WITH 2D GAUSSIAN ------------" << endl;
  if (VERBOSE)
    cout << endl;
/* Here is the definition of the Gaus2D fit parameters 
 8
 * The array gaus2dfit should be initialized by some other means 
 * before calling this function. For example by using
 * FindMoments() 
 */
/*gaus2dfit[0] = ci_;
  gaus2dfit[1] = wi_1e;
  gaus2dfit[2] = cj_;
  gaus2dfit[3] = wj_1e;
  gaus2dfit[4] = peak;
  gaus2dfit[5] = 0.1; */
  gaus2dfit_err[0] = 1e15;
  gaus2dfit_err[1] = 1e15;
  gaus2dfit_err[2] = 1e15;
  gaus2dfit_err[3] = 1e15;
  gaus2dfit_err[4] = 1e15;
  gaus2dfit_err[5] = 1e15;
  if (VERBOSE)
    cout << endl <<
      "------------ Fitting with 2D Gaussian ------------" << endl;
  if (!p->blanks)
    {
      fit2dgaus_err (columndensity, gaus2dfit, gaus2dfit_err);
      if (mott)
	{
	  /*Use regular 2d fit para as a initial guess for mott 2d gaussian fit */
	  gaus2dfit_mott[0] = gaus2dfit[0];
	  gaus2dfit_mott[1] = gaus2dfit[1];
	  gaus2dfit_mott[2] = gaus2dfit[2];
	  gaus2dfit_mott[3] = gaus2dfit[3];
	  gaus2dfit_mott[4] = gaus2dfit[4];
	  gaus2dfit_mott[5] = 1;
	  fit2dmottgaus_neldermead (columndensity, gaus2dfit_mott);
	  nfit_mott =
	    gaus2dfit_mott[4] * M_PI * gaus2dfit_mott[1] * gaus2dfit_mott[3] *
	    (1 +
	     2 / M_SQRTPI * gaus2dfit_mott[5] * exp (-gaus2dfit_mott[5] *
						     gaus2dfit_mott[5]) -
	     gsl_sf_erf (gaus2dfit_mott[5]) +
	     4.0 / 3.0 / M_SQRTPI * pow (gaus2dfit_mott[5],
					 3) * exp (-gaus2dfit_mott[5] *
						   gaus2dfit_mott[5]));
	  peakd_mott =
	    gaus2dfit_mott[4] / M_SQRTPI / pow (gaus2dfit_mott[1] *
						gaus2dfit_mott[3],
						0.5) *
	    exp (-gaus2dfit_mott[5] * gaus2dfit_mott[5]) / pow (p->magnif *
								1e-4, 3);
	  make_gaus2d_inspect (columndensity, gaus2dfit_mott,
			       p->shotnum_fileout.c_str (), true);
	}
      make_gaus2d_inspect (columndensity, gaus2dfit,
			   p->shotnum_fileout.c_str ());
      gaus2d_eval_Azimuth (gaus2dfit, p->shotnum_fileout);
    }

  if (VERBOSE)
    cout << endl;
  nfit = gaus2dfit[4] * M_PI * gaus2dfit[1] * gaus2dfit[3];
  if (VERBOSE || p->imgverbose)
    {

      printf ("\tnfit  = %.3f\n", nfit);
      printf ("\tpeakd = %.3f\n", gaus2dfit[4]);
      printf ("\tsize1 = %.3f\n", gaus2dfit[1]);
      printf ("\tsize2 = %.3f\n", gaus2dfit[3]);
    }
  nfit_err =
    pow (pow (gaus2dfit_err[4] * nfit / gaus2dfit[4], 2) +
	 pow (gaus2dfit_err[1] * nfit / gaus2dfit[1],
	      2) + pow (gaus2dfit_err[3] * nfit / gaus2dfit[3], 2), 0.5);

  //--peakd obtainied using only size along i (up to 04/16/2013)
  peakd = gaus2dfit[4] / (pow (M_PI, 0.5) * gaus2dfit[1]) / pow (p->magnif * 1e-4, 3);	// cm^-3

  //--peakd_sph  is obtained using the geometric mean
  peakd_sph = gaus2dfit[4] / (pow (M_PI * gaus2dfit[1] * gaus2dfit[3], 0.5)) / pow (p->magnif * 1e-4, 3);	// cm^-3

  if (VERBOSE)
    {
      printf ("..............  GAUSSIAN 2D FIT RESULTS ..............\n");
      printf ("ci  	  = %.1f pixels\n", gaus2dfit[0]);
      printf ("wi   	  = %.1f pixels\n", gaus2dfit[1]);
      printf ("cj    	  = %.1f pixels\n", gaus2dfit[2]);
      printf ("wj   	  = %.1f pixels\n", gaus2dfit[3]);
      printf ("peak 	  = %.3e \n", gaus2dfit[4]);
      printf ("offset	  = %.3e \n", gaus2dfit[5]);
      printf ("N from fit = %.3e \n", nfit);
      if (mott)
	{
	  printf ("wi_mott    = %.1f pixels\n", gaus2dfit_mott[1]);
	  printf ("ci_mott	  = %.1f pixels\n", gaus2dfit_mott[0]);
	  printf ("cj_mott    = %.1f pixels\n", gaus2dfit_mott[2]);
	  printf ("wj_mott    = %.1f pixels\n", gaus2dfit_mott[3]);
	  printf ("r0_mott    = %.1f pixels\n", gaus2dfit_mott[5]);
	  printf ("peak_mott  = %.3e \n", gaus2dfit_mott[4]);
	  printf ("N_mott from fit = %.3e \n", nfit_mott);
	}
    }
  TF = p->hvbar * pow (6 * nfit, 1. / 3.);
  if (VERBOSE)
    {
      printf
	("..........  DETERMINE FERMI TEMPERATURE FROM NUMBER ..........\n");
      printf (" h * ( v_r * v_r * va )^1/3 = %.3f uK\n", p->hvbar);
      printf ("                         TF = %.3f uK\n", p->hvbar);
    }

  residuals = gaus2d_residual (columndensity, gaus2dfit);
  unsigned int ci = coerce_matrix_index ((unsigned int) floor (gaus2dfit[0]),
					 columndensity->size1);
  unsigned int cj = coerce_matrix_index ((unsigned int) floor (gaus2dfit[2]),
					 columndensity->size2);
  if (VERBOSE)
    {
      cout << "ci = " << ci << " ; cj = " << cj << endl;
    }

  number_fit = 0.0;
  gsl_vector *gaus2d_v = gsl_vector_alloc (6);
  for (int e = 0; e < 6; e++)
    gsl_vector_set (gaus2d_v, e, gaus2dfit[e]);
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  //Integrate the gaus2d fit, subtracting the offset
	  number_fit += gaus2d_model (i, j, gaus2d_v) - gaus2dfit[5];
	}
    }
  return;
}

void
Fermions::FitScatt2DGauss ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT SCATTERED PHOTONS WITH 2D GAUSISAN ------------"
      << endl;
  if (VERBOSE)
    cout << endl;
  scatt2dfit[0] = gaus2dfit[0];
  scatt2dfit[1] = gaus2dfit[1];
  scatt2dfit[2] = gaus2dfit[2];
  scatt2dfit[3] = gaus2dfit[3];
  scatt2dfit[4] = maxNSP;
  scatt2dfit[5] = 0.1;
  if (VERBOSE)
    cout << endl <<
      "------------ Fitting with 2D Gaussian ------------" << endl;
  if (!p->blanks)
    fit2dgaus (columndensity_scattered_ph, scatt2dfit);
  if (VERBOSE)
    cout << endl;
  scatt2dfit[1] = fabs (scatt2dfit[1]);
  scatt2dfit[3] = fabs (scatt2dfit[3]);
  double nscattph = scatt2dfit[4] * M_PI * scatt2dfit[1] * scatt2dfit[3];
  if (VERBOSE)
    {
      printf
	("..............  SCATTERED PHOTONS 2D FIT RESULTS ..............\n");
      printf ("ci     = %.1f pixels\n", scatt2dfit[0]);
      printf ("wi     = %.1f pixels\n", scatt2dfit[1]);
      printf ("cj     = %.1f pixels\n", scatt2dfit[2]);
      printf ("wj     = %.1f pixels\n", scatt2dfit[3]);
      printf ("peak   = %.3e \n photons", scatt2dfit[4]);
      printf ("offset = %.3e \n", scatt2dfit[5]);
      printf ("N from scattered photons = %.3e \n", nscattph);	// / (p->hc * 1e3 / (p->lambda * 1e-2)) / ( ( maxI/ (1+2*maxI) ) * 2 * M_PI * 5.9e6 * p->texp ) );
    }

  return;
}

void
Fermions::FitProbe2DGauss ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT PROBE WITH 2D GAUSISAN ------------" << endl;
  if (VERBOSE)
    cout << endl;
  probe2dfit[0] = 250.;
  probe2dfit[1] = 400.;
  probe2dfit[2] = 250.;
  probe2dfit[3] = 400.;
  probe2dfit[4] = 0.75 * maxCN;
  if (VERBOSE)
    cout << endl <<
      "------------ Fitting with 2D Gaussian ------------" << endl;
  if (!p->blanks)
    fit2dgaus_no_offset (probe, probe2dfit);
  string probeinspect = p->shotnum_fileout;
  probeinspect += "_probe";
  double probefit[6] = {
    probe2dfit[0], probe2dfit[1], probe2dfit[2], probe2dfit[3],
    probe2dfit[4], 0.0
  };
  make_gaus2d_inspect (probe, probefit, probeinspect.c_str ());
  if (VERBOSE)
    cout << endl;
  probe2dfit[1] = fabs (probe2dfit[1]);
  probe2dfit[3] = fabs (probe2dfit[3]);
  if (VERBOSE)
    {
      printf ("..............  PROBE 2D FIT RESULTS ..............\n");
      printf ("ci     = %.1f pixels\n", probe2dfit[0]);
      printf ("wi     = %.1f pixels\n", probe2dfit[1]);
      printf ("cj     = %.1f pixels\n", probe2dfit[2]);
      printf ("wj     = %.1f pixels\n", probe2dfit[3]);
      printf ("peak   = %.3e counts\n", probe2dfit[4]);
    }

  printf (" 1/e^2 beam waist of the probe:\n");
  printf ("  i (horizontal on camera) = %.3f um\n", probe2dfit[1] * sqrt (2) * 16.);	//sqrt(2) to go from 1/e to 1/e^2 and 16um/px for the Andor
  printf ("  j (  vertical on camera) = %.3f um\n", probe2dfit[3] * sqrt (2) * 16.);	//sqrt(2) to go from 1/e to 1/e^2 and 16um/px for the Andor
  unsigned int smooth_bins = 10;
  gsl_matrix *probe_smoothed = smooth (probe, smooth_bins);
  unsigned int probemax_i, probemax_j;
  double probe_max_pos;
  findpeak (probe_smoothed, &probemax_i, &probemax_j, &probe_max_pos, true);
  printf ("        max probe smoothed = %.3f counts\n\n", probe_max_pos);
  printf ("             max probe fit = %.3f counts\n", probe2dfit[4]);
  return;
}


void
Fermions::Fit2DFermi ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT COLUMN DENSITY WITH 2D FERMI-DIRAC ------------"
      << endl;
  fermi2dfit[0] = gaus2dfit[4];
  fermi2dfit[1] = -5.0;
  fermi2dfit[2] = gaus2dfit[1];
  fermi2dfit[3] = gaus2dfit[3];
  fermi2dfit[4] = gaus2dfit[0];
  fermi2dfit[5] = gaus2dfit[2];
  fermi2dfit[6] = gaus2dfit[5];
  if (VERBOSE)
    {
      const char *const pnames[] = {
	"peakd", "BetaMu", "radius_i", "radius_j", "center_i",
	"center_j", "offset"
      };
      printf ("Starting fit values:\n");
      for (int e = 0; e < 7; e++)
	{
	  printf (" fermi2dfit[%d] = %f\t%s\n", e, fermi2dfit[e], pnames[e]);
	}
    }

  //override so that it doesn't take forever
/*  fermi2dfit[0] = 1.63e+02;
  fermi2dfit[1] = 1.76e+00;
  fermi2dfit[2] = 68.4;
  fermi2dfit[3] = 87.9;
  fermi2dfit[4] = 265.0;
  fermi2dfit[5] = 297.4;
  fermi2dfit[7] = 1.38; */
  if (!p->blanks)
    {
      fit2dfermi_neldermead (columndensity, fermi2dfit);
      make_fermi2d_gaus2d_inspect (columndensity, fermi2dfit,
				   gaus2dfit, p->shotnum_fileout.c_str ());
      fermi2d_eval_Azimuth (fermi2dfit, p->shotnum_fileout);
    }

  if (VERBOSE)
    cout << endl;
  BetaMu = fermi2dfit[1];
  TF_2d = pow (6 * f2 (BetaMu), -0.333);
  Fugacity_Fermi = exp (BetaMu);
  f_Fermi = fq (BetaMu);
  ri_Fermi = fermi2dfit[2];
  rj_Fermi = fermi2dfit[3];
  T_2d_rd =
    p->B1 / (2 * p->kbm) * ri_Fermi * ri_Fermi * p->magnif *
    p->magnif / f_Fermi;
  T_2d_ax =
    p->B2 / (2 * p->kbm) * rj_Fermi * rj_Fermi * p->magnif *
    p->magnif / f_Fermi;
  if (VERBOSE || p->showfermi)
    {
      printf ("..............  FERMI 2D FIT RESULTS ..............\n");
      printf ("n0     = %.2e\n", fermi2dfit[0]);
      printf
	("BetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	 fermi2dfit[1], f2 (fermi2dfit[1]), TF_2d);
      printf ("ri     = %.1f pixels\n", fermi2dfit[2]);
      printf ("rj     = %.1f pixels\n", fermi2dfit[3]);
      printf ("ci     = %.1f pixels\n", fermi2dfit[4]);
      printf ("cj     = %.1f pixels\n", fermi2dfit[5]);
      printf ("B      = %.2e \n", fermi2dfit[6]);
      printf ("fugacity   = %.2e\n", Fugacity_Fermi);
      printf ("f(z)       = %.2e\n", f_Fermi);
      printf ("T_2d_rd   = %.2f uK\n", T_2d_rd);
      printf ("T_2d_ax   = %.2f uK\n", T_2d_ax);
      printf ("TF     = %.2f uK\n", TF);
    }

  n0 = fermi2dfit[0];
  ci_Fermi = fermi2dfit[4];
  cj_Fermi = fermi2dfit[5];
  B_Fermi = fermi2dfit[6];
  return;
}



void
Fermions::GetAzimuthalAverageEllipse ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ GET AZIMUTHAL AVERAGE OF CLOUD (ELLIPSE) ------------"
      << endl;
  //Aspect ratio from trap geometry and self similar expansion
  double trap_aspect = sqrt (p->B2 / p->B1);
  //Aspect ratio obtained from 2D Gaussian fit
  double gaus2d_aspect = gaus2dfit[3] / gaus2dfit[1];
  // Aspect ratio determined from trap geometry and self similar expansion
  // should not be different by more than 2%
  // If they are it means that we do not know the trap frequencies
  double ar_err = abs (trap_aspect - gaus2d_aspect) / trap_aspect;

  //Disabled aspect ratio warning message, since we now have a
  //variety of traps that we image and they have different AR's
  if (ar_err > 0.02 && false)
    {
      printf
	("\n******   WARNING: 2D-Gauss aspect ratio != Self similar aspect ratio !!! ******\n");
      cout << "\tAspect ratio from trap geometry = " << trap_aspect << endl;
      cout << "\tAspect ratio from 2D Gauss = " << gaus2d_aspect << endl;
      cout << "\tDiscrepancy = " << ar_err * 100 << " %" << endl;
    }

  // Define aspect ratio to be used in calculating the azimuthal average
  p->AR = gaus2d_aspect;
  if (VERBOSE)
    {
      cout << "\tAspect ratio from trap geometry = " << trap_aspect << endl;
      cout << "\tAspect ratio from 2D Gauss = " << gaus2d_aspect << endl;
      cout << "\tAspect ratio used = " << p->AR << endl;
    }

  // Define the number of bins to be used in the radial profile  and the bin size 
  nbins = 2048;
  binsize = 1.2;		//in pixels
  // Declare the auxiliary arrays to calculate the average 
  //    _0 is sum of values for the data
  //    _1 is number of values 
  //    average at bin =  _0 / _1
  //
  gsl_vector *azim_histo_0 = gsl_vector_alloc (nbins);
  gsl_vector *azim_histo_1 = gsl_vector_alloc (nbins);
  // Define a matrix that stores the reduced radius rho for each pixel
  gsl_matrix *rho = gsl_matrix_alloc (columndensity->size1,
				      columndensity->size2);
  // Initialize arrays
  for (unsigned int j = 0; j < nbins; j++)
    {
      gsl_vector_set (azim_histo_0, j, 0.0);
      gsl_vector_set (azim_histo_1, j, 0.0);
    }

  if (VERBOSE)
    {
      printf ("\t...Begin iteration for azimuthal average calculation\n");
    }
  // initialize cut arrays
  icut_r = gsl_vector_alloc (columndensity->size1);
  icut_dat = gsl_vector_alloc (columndensity->size1);
  jcut_r = gsl_vector_alloc (columndensity->size2);
  jcut_dat = gsl_vector_alloc (columndensity->size2);
  bool DEBUG_VECTORS = false;
  if (DEBUG_VECTORS)
    {
      printf ("icut size = %u\n", (unsigned int) columndensity->size1);
      printf ("jcut size = %u\n", (unsigned int) columndensity->size2);
    }

  bool DEBUG_AZIMUTH = false;
  // populate auxiliary arrays for average
  // at the same time populate arrays for cuts along principal axis 
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  double xi = (double) i;	// this is y  radial
	  double xj = (double) j;	// this is g  axial
	  if (DEBUG_AZIMUTH)
	    printf ("i=%u,j=%u\n", i, j);
	  double data = gsl_matrix_get (columndensity, i, j);
	  if (DEBUG_AZIMUTH)
	    printf ("data = %f\n", data);
	  // Find the radial distance, aspect ratio is taken into account
	  unsigned int dist =
	    (unsigned int) floor (pow
				  (pow (p->AR * (xi - gaus2dfit[0]), 2) +
				   pow (xj - gaus2dfit[2], 2),
				   0.5) / binsize);
	  if (dist > nbins)
	    {
	      printf ("\n!!!!!!!!!!!!!!!!!!!\n");
	      printf
		(" ERROR in AZIMUTHAL AVERAGE CALCULATION:\n\t dist > nbins\n");
	      printf ("!!!!!!!!!!!!!!!!!!!\n");
	      exit (0);
	    }

	  gsl_matrix_set (rho, i, j, dist);
	  if (DEBUG_AZIMUTH)
	    printf ("rho = %u\n", dist);
	  if (DEBUG_AZIMUTH)
	    printf (" dist(%d,%d) from (%.2f,%.2f)  = %d\n", i, j,
		    gaus2dfit[0], gaus2dfit[2], dist);
	  //increment sum of data
	  gsl_vector_set (azim_histo_0, dist,
			  gsl_vector_get (azim_histo_0, dist) + data);
	  //increment N
	  gsl_vector_set (azim_histo_1, dist,
			  gsl_vector_get (azim_histo_1, dist) + 1);
	  //set icut
	  if (j == (unsigned int) floor (gaus2dfit[2]))
	    {
	      gsl_vector_set (icut_r, i, p->AR * (i - gaus2dfit[0]));
	      gsl_vector_set (icut_dat, i, data);
	    }

	  //set jcut
	  if (i == (unsigned int) floor (gaus2dfit[0]))
	    {
	      gsl_vector_set (jcut_r, j, j - gaus2dfit[2]);
	      gsl_vector_set (jcut_dat, j, data);
	    }


	}
    }

  char base[MAXPATHLEN];
  getcwd (base, MAXPATHLEN);
  string rho_path = makepath (base, p->shotnum_fileout, "_rho.TIFF");
  string rho_ascii_path = makepath (base, p->shotnum_fileout, "_rho.ascii");
  if (p->savetiff)
    toTiffImage (rho, rho_path, true);
  if (p->saveascii)
    save_gsl_matrix_ASCII (rho, rho_ascii_path);
  if (p->azimverbose)		//CAN CHECK HERE THE INTEGRITY OF THE AZIMUTHAL AVERAGE HISTOGRAMS
    {
      printf ("\n AZIMUTHAL AVERAGE HISTOGRAMS:\n");
      for (unsigned int ii = 0; ii < nbins; ii++)
	{
	  printf (" azim_histo_0( i = %d ) = %.3f ", ii,
		  gsl_vector_get (azim_histo_0, ii));
	  printf (" azim_histo_1( i = %d ) = %.3f \n", ii,
		  gsl_vector_get (azim_histo_1, ii));
	}
    }


  // Determine how many bins were used on the auxiliary arrays
  usedbins = nbins - 1;
  while (gsl_vector_get (azim_histo_1, usedbins) == 0)
    usedbins--;
  // Fill out arrays that contain all the data
  azimuthal_all_r = gsl_vector_alloc (usedbins);
  azimuthal_all_dat = gsl_vector_alloc (usedbins);
  double r;
  if (VERBOSE)
    printf ("\n AZIMUTHAL AVERAGE RADIUS AND DATA VECTORS:\n");
  for (unsigned int index = 0; index < usedbins; index++)
    {
      r = index * binsize;
      double hsum = gsl_vector_get (azim_histo_0, index);
      double hnum = gsl_vector_get (azim_histo_1, index);
      if (VERBOSE)
	{
	  if (hnum == 0)
	    {
	      printf
		("ERROR: no entries in azimuthal histogram at radius r = %.6f\n",
		 r);
	    }
	  if (false)
	    {
	      printf (" r = %.3f,  az_avg(r) = %.3f\n", r, hsum / hnum);
	    }

	}
      gsl_vector_set (azimuthal_all_r, index, r);
      gsl_vector_set (azimuthal_all_dat, index, hsum / hnum);
    }


  // Fill out arrays that contain cropped data for the fit
  fitbins = usedbins;
  // Azimuthal average should not extend to a radius larger than azimuth_maxr
  if (usedbins * binsize > p->azimuth_maxr)
    {
      fitbins = (unsigned int) ceil (p->azimuth_maxr / binsize);
    }

  // Also one can chop the trailing values of the average
  else
    {
      fitbins -= (unsigned int) ceil (p->azimuth_chop / binsize);
    }


  // And also one can take away values near thecenter
  unsigned int start = 0;
  if (p->azimuth_start > 0.)
    {
      start = (unsigned int) ceil (p->azimuth_start / binsize);
    }

  // After computing the sums, divide by the number to get the azimuthal average
  // 0 is the sum of data
  // 1 is the number of pixels included in the sum

  azimuthal_r = gsl_vector_alloc (fitbins - start);
  azimuthal_dat = gsl_vector_alloc (fitbins - start);
  for (unsigned int index = start; index < fitbins; index++)
    {
      gsl_vector_set (azimuthal_r, index - start,
		      gsl_vector_get (azimuthal_all_r, index));
      gsl_vector_set (azimuthal_dat, index - start,
		      gsl_vector_get (azimuthal_all_dat, index));
    }


  to_dat_file_2 (azimuthal_r, azimuthal_dat, p->shotnum_fileout,
		 "datAzimuth.AZASCII");
  to_dat_file_2 (azimuthal_all_r, azimuthal_all_dat,
		 p->shotnum_fileout, "datAllAzimuth.AZASCII");
  to_dat_file_2 (icut_r, icut_dat, p->shotnum_fileout, "datIcut.AZASCII");
  to_dat_file_2 (jcut_r, jcut_dat, p->shotnum_fileout, "datJcut.AZASCII");
  return;
  /* Azimuthal-type data consists of two vector with an equal number of 
   * elements.  Each element in the vectors represents a bin. 
   * One vector contains the distances from the center and the other
   * vector contains the averaged value.  
   *
   */
}


void
Fermions::FitAzimuthalFermi ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT AZIMUTHAL AVERAGE WITH FERMI-DIRAC ------------"
      << endl;
  /********************************* Finite temperature azimuthal fit ***/
  fermi_azimuth_fit[0] = gaus2dfit[4];
  fermi_azimuth_fit[1] = -5.0;
  //fermi_azimuth_fit[2] = pow(wi_1e*wj_1e,0.5); 
  fermi_azimuth_fit[2] = gaus2dfit[3];
  fermi_azimuth_fit[3] = gaus2dfit[5];
  fermi_azimuth_fit[4] = 0.1;
  if (VERBOSE)
    {
      const char *const pnames[] = {
	"peakd", "BetaMu", "radius", "offset", "slope"
      };
      printf ("Starting fit values:\n");
      for (int e = 0; e < 5; e++)
	{
	  printf (" fermi_azimuth_fit[%d] = %f\t%s\n", e,
		  fermi_azimuth_fit[e], pnames[e]);
	}
    }

  gsl_vector *azimuthal_[2] = {
    azimuthal_r, azimuthal_dat
  };
  if (!p->blanks)
    {
      fit1dfermi_azimuthal_neldermead (azimuthal_, fermi_azimuth_fit);
      if (VERBOSE || p->showfermi)
	{
	  printf
	    ("..............  FERMI AZIMUTHAL FIT RESULTS ..............\n");
	  printf ("n0         = %.2e\n", fermi_azimuth_fit[0]);
	  printf ("BetaMu     = %.2e\n", fermi_azimuth_fit[1]);
	  printf ("r          = %.1f pixels\n", fermi_azimuth_fit[2]);
	  printf ("offset     = %.1f pixels\n", fermi_azimuth_fit[3]);
	  printf ("mx         = %.1f pixels\n", fermi_azimuth_fit[4]);
	}
      gsl_vector *fitAzimuth =
	fermiAzimuth_eval (azimuthal_r, fermi_azimuth_fit);
      to_dat_file_2 (azimuthal_r, fitAzimuth, p->shotnum_fileout,
		     "fitAzimuth.AZASCII");
    }

  fermi_azimuth_fit[2] = fabs (fermi_azimuth_fit[2]);
  r_az = fermi_azimuth_fit[2];
  TF_az = pow (6 * f2 (fermi_azimuth_fit[1]), -0.333);
  BetaMu_az = fermi_azimuth_fit[1];
  Fugacity_az = exp (BetaMu_az);
  f_az = fq (BetaMu_az);
  T_az = p->B1 / (2 * p->kbm) * r_az * r_az * p->magnif * p->magnif / f_az;
  n0_az = fermi_azimuth_fit[0];
  B_az = fermi_azimuth_fit[3];
  mx_az = fermi_azimuth_fit[4];
  if (VERBOSE)
    cout << endl <<
      "---------- FIT AZIMUTHAL AVERAGE WITH T=0 FERMI-DIRAC -----------"
      << endl;
  /********************************* Zero temperature azimuthal fit *****/
  fermi_azimuth_fit_zero[0] = n0_az;
  fermi_azimuth_fit_zero[1] = r_az;
  fermi_azimuth_fit_zero[2] = B_az;
  fermi_azimuth_fit_zero[3] = mx_az;
  if (!p->blanks)
    {
      fit1dfermi_azimuthal_zero_neldermead (azimuthal_,
					    fermi_azimuth_fit_zero);
      if (VERBOSE || p->showfermi)
	{
	  printf
	    ("............  FERMI AZIMUTHAL T=0 FIT RESULTS .............\n");
	  printf ("n0         = %.2e\n", fermi_azimuth_fit_zero[0]);
	  printf ("BetaMu     = %.2e\n", fermi_azimuth_fit_zero[1]);
	  printf ("r          = %.1f pixels\n", fermi_azimuth_fit_zero[2]);
	  printf ("offset     = %.1f pixels\n", fermi_azimuth_fit_zero[3]);
	  printf ("mx         = %.1f pixels\n", fermi_azimuth_fit_zero[4]);
	}
      gsl_vector *fitAzimuthZeroT = fermiAzimuthZeroT_eval (azimuthal_r,
							    fermi_azimuth_fit_zero);
      to_dat_file_2 (azimuthal_r, fitAzimuthZeroT, p->shotnum_fileout,
		     "fitAzimuthZeroT.AZASCII");
    }

  fermi_azimuth_fit[1] = fabs (fermi_azimuth_fit[1]);
  n0_az_zeroT = fermi_azimuth_fit_zero[0];
  r_az_zeroT = fermi_azimuth_fit_zero[1];
  B_az_zeroT = fermi_azimuth_fit_zero[2];
  mx_az_zeroT = fermi_azimuth_fit_zero[3];
  if (VERBOSE)
    cout << endl;
  if (VERBOSE || p->showfermi)
    {
      printf ("..............  FERMI AZIMUTHAL FIT RESULTS ..............\n");
      printf ("n0         = %.2e\n", fermi_azimuth_fit[0]);
      printf
	("BetaMu     = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	 fermi_azimuth_fit[1], f2 (fermi_azimuth_fit[1]), TF_az);
      printf ("r          = %.1f pixels\n", fermi_azimuth_fit[2]);
      printf ("r (zero T) = %.1f pixels\n", fermi_azimuth_fit_zero[1]);
      printf ("offset     = %.1f pixels\n", fermi_azimuth_fit[3]);
      printf ("mx         = %.1f pixels\n", fermi_azimuth_fit[4]);
      printf ("fugacity   = %.2e\n", Fugacity_az);
      printf ("f(z) = %.2e\n", f_az);
      printf ("T_az = %.2f uK\n", T_az);
      printf ("TF   = %.2f uK\n", TF);
    }



  return;
}


void
Fermions::MakePlots ()
{
  if (VERBOSE)
    cout << endl << "----------------- MAKE PLOTS ----------------" << endl;
  stringstream inspectstr;
  inspectstr << "inspectAz_ascii.py ";
  inspectstr << p->shotnum;
  if (p->andor2)
    inspectstr << " andor2";
  system (inspectstr.str ().c_str ());
  return;
}


void
Fermions::ComputeIntegrated1DDensity ()
{

  if (VERBOSE)
    cout << endl <<
      "----------- COMPUTE INTEGRATED 1D DENSITY ------------" << endl;
  unsigned int s1 = columndensity->size1;
  unsigned int s2 = columndensity->size2;
  sum_density_0 = gsl_vector_alloc (s1);
  sum_density_1 = gsl_vector_alloc (s2);
  sum_missing_0 = gsl_vector_alloc (s1);
  sum_missing_1 = gsl_vector_alloc (s2);
  gsl_vector_set_all (sum_density_0, 0.0);
  gsl_vector_set_all (sum_density_1, 0.0);
  gsl_vector_set_all (sum_missing_0, 0.0);
  gsl_vector_set_all (sum_missing_1, 0.0);
  double cd_ij = 0.0, mi_ij;
  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  cd_ij = gsl_matrix_get (columndensity, i, j);
	  mi_ij = gsl_matrix_get (missing_counts, i, j);
	  gsl_vector_set (sum_density_0, i,
			  gsl_vector_get (sum_density_0, i) + cd_ij);
	  gsl_vector_set (sum_density_1, j,
			  gsl_vector_get (sum_density_1, j) + cd_ij);
	  gsl_vector_set (sum_missing_0, i,
			  gsl_vector_get (sum_missing_0, i) + mi_ij);
	  gsl_vector_set (sum_missing_1, j,
			  gsl_vector_get (sum_missing_1, j) + mi_ij);
	}
    }

  if (VERBOSE)
    cout << endl;
  if (VERBOSE)
    cout << endl <<
      "------------- FIT INTEGRATED 1D WITH GAUSSIAN  --------------" << endl;
  fit1d_gaus_0[0] = gaus2dfit[0];
  fit1d_gaus_0[1] = gaus2dfit[1];
  fit1d_gaus_0[2] = gsl_vector_max (sum_density_0);
  fit1d_gaus_0[3] = 0.1;
  if (!p->blanks)
    fit1dgaus (sum_density_0, fit1d_gaus_0);
  fit1d_gaus_1[0] = gaus2dfit[2];
  fit1d_gaus_1[1] = gaus2dfit[3];
  fit1d_gaus_1[2] = gsl_vector_max (sum_density_1);
  fit1d_gaus_1[3] = 0.1;
  if (!p->blanks)
    fit1dgaus (sum_density_1, fit1d_gaus_1);
  if (VERBOSE)
    cout << endl <<
      "----------- FIT INTEGRATED 1D WITH FERMI-DIRAC  ------------" << endl;
// 02/07/2012 Commented out fermi fits on axial and radial profiles
  fit1d_fermi_0[0] = gsl_vector_max (sum_density_0);
  fit1d_fermi_0[1] = -5.0;
  fit1d_fermi_0[2] = gaus2dfit[1];
  fit1d_fermi_0[3] = gaus2dfit[0];
  fit1d_fermi_0[4] = fit1d_gaus_0[3];
  if (p->fitfermi1D && !p->blanks)
    {
      fit1dfermi_neldermead (sum_density_0, fit1d_fermi_0);
    }

  if (VERBOSE)
    printf ("\n---> Finished _0 axis\n\n");
  fit1d_fermi_1[0] = gsl_vector_max (sum_density_1);
  fit1d_fermi_1[1] = -5.0;
  fit1d_fermi_1[2] = gaus2dfit[3];
  fit1d_fermi_1[3] = gaus2dfit[2];
  fit1d_fermi_1[4] = fit1d_gaus_1[3];
  if (p->fitfermi1D && !p->blanks)
    fit1dfermi_neldermead (sum_density_1, fit1d_fermi_1);
  if (VERBOSE && p->fitfermi1D)
    cout << endl << "---> Calculating T/TF from fit results" << endl;
  TF_rd = pow (6 * f2 (fit1d_fermi_0[1]), -0.333);
  TF_ax = pow (6 * f2 (fit1d_fermi_1[1]), -0.333);
  if (VERBOSE)
    {
      cout << endl <<
	"------------ INTEGRATED 1D GAUS FIT RESULTS ------------" << endl;
      printf ("\n_0 Profile:\n");
      printf ("c_ref = (%.0f,%.0f)\n", abs_ci, abs_cj);
      printf ("\tc_rd = %.1f\n", fit1d_gaus_0[0]);
      printf ("\tw_rd = %.2f\n", fit1d_gaus_0[1]);
      printf ("\tA_rd = %.3e\n", fit1d_gaus_0[2]);
      printf ("\tB_rd = %.3e\n", fit1d_gaus_0[3]);
      printf ("\n_1 Profile:\n");
      printf ("\tc_ax = %.1f\n", fit1d_gaus_1[0]);
      printf ("\tw_ax = %.2f\n", fit1d_gaus_1[1]);
      printf ("\tA_ax = %.3e\n", fit1d_gaus_1[2]);
      printf ("\tB_ax = %.3e\n", fit1d_gaus_1[3]);
    }

  if (VERBOSE && p->fitfermi1D)
    {

      cout << endl <<
	"------------ INTEGRATED 1D DENSITY FERMI FIT RESULTS ------------"
	<< endl;
      printf ("\n_0 Profile:\n");
      printf ("\tn0_rd     = %.1f\n", fit1d_fermi_0[0]);
      printf
	("\tBetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	 fit1d_fermi_0[1], f2 (fit1d_fermi_0[1]), TF_rd);
      printf ("\tr_rd      = %.3e\n", fit1d_fermi_0[2]);
      printf ("\tc_rd      = %.3e\n", fit1d_fermi_0[3]);
      printf ("\tB_rd      = %.3e\n", fit1d_fermi_0[4]);
      printf ("\n_1 Profile:\n");
      printf ("\tn0_ax     = %.1f\n", fit1d_fermi_1[0]);
      printf
	("\tBetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	 fit1d_fermi_1[1], f2 (fit1d_fermi_1[1]), TF_ax);
      printf ("\tr_ax      = %.3e\n", fit1d_fermi_1[2]);
      printf ("\tc_ax      = %.3e\n", fit1d_fermi_1[3]);
      printf ("\tB_ax      = %.3e\n", fit1d_fermi_1[4]);
    }


  return;
}
