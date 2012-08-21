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

extern bool VERBOSE;


struct params
{
  unsigned int shot;

  string shotnum, reportfile, atomsfile, noatomsfile, atomsreffile,
    noatomsreffile;

  // ROI is defined as (ax0pos, ax1pos, ax0size, ax1size)
  unsigned int roi[4], roisize[2], centerpt[2];
  bool keeproi;

  bool verbose, center, crop, plots, reanalyze, roi_user,
    fitfermi1D, fermi2d, fermiazimuth, showfermi, show_B, blanks;

  // Whether or not to save ASCII or TIFF files
  // Useful for debugging the algorithms
  bool saveascii, savetiff;


  double azimuth_maxr, azimuth_chop, azimuth_start;

  double lambda, hc, gamma, magnif, kbm, decay;

  double h;

  // Image parameters
  double texp, odttof, det;

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
  p->det = (double) getINI_num (p->reportfile, "ANDOR", "phcdet");	// phase contrast detuning in MHz
  p->phc = (bool) getINI_num (p->reportfile, "ANDOR", "phc");	// phase contrast 
  p->phcsign = (double) getINI_num (p->reportfile, "ANDOR", "phcsign");	// prefactor in the calculation of the phase contrast signal

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


  p->lambda = 670.977e-7;	// cm
  p->hc = 1.98644521e-25;	// J*m 
  p->gamma = 5.9e6;		// Hertz
  p->decay = 27.2e-9;		// Seconds

  p->mJ_per_photon = 1.e3 * (p->hc / (p->lambda * 1.e-2));
  p->photon_per_mJ = 1. / p->mJ_per_photon;


  p->kbm = 1385.;		// um^2/ms^2/uK This is kb/m

  double emp = 0.88;
  p->andoreff_10MHz_14bit_x1_Electron_Mult =
    p->mJ_per_photon * (67.9 / 0.95 / emp);

  // mJ/count = (mJ per photon) * (electrons per A/D count  / QuantumEff / EmpiricalCorrectionFactor) 
  emp = 0.95;
  p->andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset =
    p->mJ_per_photon * (12.9 / 0.95 / emp);

  //p->eff = p->andoreff_10MHz_14bit_x1_Electron_Mult;
  p->eff = p->andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset;

  p->magnif = 16. / 5.;		// 16um/pixel for the Andor, using the 5x obj with the telephoto at f=200


  p->h = 48.;			// uK/MHz  this is Planck's constant divided by Boltzmann's constant 
  // Obtain trap depth from report
  p->trapdepth = p->finalcpow / 10. * p->odtmaxdepth;
  // Calculate the geometric mean of trap frequencies
  // Useful to obtain the Fermi temperature
  // EF = h vbar (6N)^(1/3)
  p->hvbar = p->h * pow (p->odtv0radial * p->odtv0radial * p->odtv0axial, 1. / 3.) * sqrt (p->finalcpow / 10.) * 1e-6;	// 1e-6 is to give the trap freqs in MHz
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
  void Fit2DGauss ();
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
  double nfit_err;
  double peakd;

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

  atoms = ReadFitsImg_gsl_matrix (p->atomsfile);
  noatoms = ReadFitsImg_gsl_matrix (p->noatomsfile);
  atomsref = ReadFitsImg_gsl_matrix (p->atomsreffile);
  noatomsref = ReadFitsImg_gsl_matrix (p->noatomsreffile);


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
  string norm_path = makepath (base, p->shotnum, "_normpatch.TIFF");
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
  eigenfacestr << pos0 << "," << pos1 << "," << siz0 << "," << siz1;
  //printf ("%s\n", eigenfacestr.str ().c_str ());
  //system (eigenfacestr.str ().c_str ());


  return;
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
  double maxT1, maxT2;		//maxima for term1 and term2 in the column density
  maxT1 = -1e6;
  maxT2 = -1e6;

  double sigma0 = 3 * M_PI * pow (p->lambda / 2 / M_PI, 2);	//This corresponds to half the maximal cross-section, because of our polarization.  
  double isat = 10.2;		// This corresponds to twice the minimal Isat because of our polarization

  //should be equal to  gamma h v / sigma0 .... test this.  


  double eff = p->eff;

  double det = 0.;		//for now all  absorption imaging is on resonance 
  maxI = 0.;

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



  /************ PHASE-CONTRAST IMAGING VARIABLES ************/
  double signal, phase_shift, sig1, sig2;
  gsl_matrix *phc_signal;
  if (p->phc)
    {
      phc_signal = gsl_matrix_alloc (s1, s2);	// A matrix of the phase-contrast Signal at each pixel.
      sigma0 = 6 * M_PI * pow (p->lambda / 2. / M_PI, 2);	//This is the maximal cross section
      isat = 5.1;		// This is the standard saturation intensity
      det = p->det * 1.e6 / p->gamma;	//detuning in units of gamma
      cout << endl;
      cout << "   Andor magnification     = " << p->magnif << " um/pixel" <<
	endl;
      cout << "   Phase-contrast detuning = " << det << " Gamma  = " << det *
	p->gamma / 1.e6 << " MHz" << endl;
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
	  i0 = (c0 * eff / (pow (p->magnif * 1e-4, 2)) / p->texp / (isat));	// intensity in noatoms frame in units of isat
	  i0eff = (c0 * eff / (pow (p->magnif * 1e-4, 2)) / p->texp / (isat * p->alphastar));	// effective intensity in noatoms frame
	  if (i0 > maxI)
	    maxI = i0;

	  OD = log (fabs (c0 / c1));

	  if (OD > OD_err)
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

	  if (!p->phc)
	    {

	      term1 = p->alphastar * (1. + 4. * det * det) * OD;


	      pee = i0eff / (1 + 2 * i0eff);
	      if (pee > maxPEE)
		maxPEE = pee;
	      //sp = scattered photons in pixel
	      sp = (c0 - c1) * eff * p->photon_per_mJ;
	      if (sp > maxSP)
		maxSP = sp;
	      //nsp = number of atoms in pixel, calculated from number of scattered photons
	      nsp = sp / (pee * p->texp / p->decay);
	      if (nsp > maxNSP)
		maxNSP = nsp;

	      Tsp += sp;	// total number of scattered photons
	      Nsp += nsp;	// total number of atoms from scattered photons 

	      term2 = 2 * i0 * (1. - c1 / c0);

	      if (term1 > maxT1)
		maxT1 = term1;
	      if (term2 > maxT2)
		maxT2 = term2;

	      cd = (term1 + term2) / sigma0 * pow (p->magnif * 1e-4, 2);

	    }

	  /*********** PHASE - CONTRAST IMAGING ***************/

	  else
	    {
	      signal = p->phcsign * (c1 / c0 - 1);
	      gsl_matrix_set (phc_signal, i, j, signal);
	      phase_shift = -1.0 * det * signal / (0.5 - det);

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

	      i0 = 2. * i0;	// The polarizer cuts it down by a factor of 2.  

	      sig1 =
		(0.5 - det) * sigma0 / pow (p->magnif * 1e-4,
					    2) / (1 + 4 * det * det + 2 * i0);
	      sig2 =
		(0.5 - (det + 13.3)) * sigma0 / pow (p->magnif * 1e-4,
						     2) / (1 + 4 * (det +
								    13.3) *
							   (det + 13.3) +
							   2 * i0);

	      cd = 2 * signal / (sig1 + sig2);	//columdensity 
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
	    maxCD = cd;

	  gsl_matrix_set (columndensity_scattered_ph, i, j, nsp);
	  gsl_matrix_set (columndensity, i, j, cd);
	  gsl_matrix_set (missing_counts, i, j, c1 - c0);
	  gsl_matrix_set (probe, i, j, c0);
	}
    }

  maxT1 = maxT1 / sigma0 * pow (p->magnif * 1e-4, 2);
  maxT2 = maxT2 / sigma0 * pow (p->magnif * 1e-4, 2);

  if (p->phc)
    {
      unsigned int smooth_bins = 4;
      gsl_matrix *phc_signal_smoothed = smooth (phc_signal, smooth_bins);
      unsigned int phcsig_i, phcsig_j;
      double phcsig_max_pos, phcsig_max_neg;
      findpeak (phc_signal_smoothed, &phcsig_i, &phcsig_j, &phcsig_max_pos,
		true);
      findpeak (phc_signal_smoothed, &phcsig_i, &phcsig_j, &phcsig_max_neg,
		false);
      maxPHCSIG =
	fabs (phcsig_max_pos) >
	fabs (phcsig_max_neg) ? phcsig_max_pos : phcsig_max_neg;
    }



  if (VERBOSE)
    {
      cout << endl << "------------ Column Density Stats ------------" <<
	endl;
      if (!p->phc)
	{
	  cout << endl << "----> Method used : absorption" << endl;
	  cout << "\ts1=" << s1 << ", s2=" << s2 << endl;
	  cout << "\talphastar = " << p->alphastar << endl;
	  cout << "\tmax term1 = " << maxT1 << endl;
	  cout << "\tmax term2 = " << maxT2 << endl;
	  cout << "\tmax OD = " << maxOD << endl;
	  cout << "\tmax probe intensity = " << maxI << " Isat " << endl;
	  cout << "\tmax CD = " << maxCD << endl;
	  cout << "\tmax scattered photons = " << maxSP << endl;
	  cout << "\tmax number from scattered photons = " << maxNSP << endl;
	  cout << "\tmax rho_{ee} = " << maxPEE << endl;
	  printf ("\tcounts in atoms picture:     ( %f to %f )\n", minCA,
		  maxCA);
	  printf ("\tcounts in no atoms picture:  ( %f to %f )\n", minCN,
		  maxCN);
	  cout << "\ttotal scattered photons = " << Tsp << endl;
	  cout << "\ttotal number from scattered photons = " << Nsp << endl;
	  cout << "\ntotal number of photons in probe pulse = " << Tp0 <<
	    endl;
	  cout << "First few elements of column density matrix: " << endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (columndensity, i, j));
		}
	      cout << "... " << endl;
	    }

	}
      else
	{
	  cout << endl << "----> Method used : phase contrast" << endl;
	  cout << "\ts1=" << s1 << ", s2=" << s2 << endl;
	  cout << "\tmax probe intensity = " << maxI << " Isat " << endl;
	  cout << "\tmax PHI = " << maxPHI << endl;
	  cout << "\tmax CD = " << maxCD << endl;
	  printf ("\tcounts in atoms picture:     ( %f to %f )\n", minCA,
		  maxCA);
	  printf ("\tcounts in no atoms picture:  ( %f to %f )\n", minCN,
		  maxCN);
	  cout << "\ntotal number of photons in probe pulse = " << Tp0 <<
	    endl;
	  cout << "First few elements of column density matrix: " << endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (columndensity, i, j));
		}
	      cout << "... " << endl;
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

  string column_path = makepath (base, p->shotnum, "_column.TIFF");
  string scatt_path = makepath (base, p->shotnum, "_column_scatt.TIFF");
  string missing_path = makepath (base, p->shotnum, "_missing_counts.TIFF");
  string probe_path = makepath (base, p->shotnum, "_probe.TIFF");

  string column_ascii_path = makepath (base, p->shotnum, "_column.ascii");
  string scatt_ascii_path =
    makepath (base, p->shotnum, "_column_scatt.ascii");
  string missing_ascii_path =
    makepath (base, p->shotnum, "_missing_counts.ascii");
  string probe_ascii_path = makepath (base, p->shotnum, "_probe.ascii");

  if (p->savetiff)
    {
      toTiffImage (columndensity, column_path);
      toTiffImage (columndensity_scattered_ph, scatt_path);
      toTiffImage (probe, probe_path);
      toTiffImage (missing_counts, missing_path);
    }

  if (p->saveascii)
    {
      save_gsl_matrix_ASCII (columndensity_scattered_ph, scatt_ascii_path);
      save_gsl_matrix_ASCII (columndensity, column_ascii_path);
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
  Gaus2DGuess (columndensity, gaus2dfit, p->shotnum, false);
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
      "------------ CROPPING COLUMN DENSITY TO MINIMAL AREA------------" <<
      endl;

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
      printf ("    Determined ROI for minimal crop [%d,%d,%d,%d]\n", roi[0],
	      roi[1], roi[2], roi[3]);
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
Fermions::Fit2DGauss ()
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
    cout << endl << "------------ Fitting with 2D Gaussian ------------" <<
      endl;

  if (!p->blanks)
    {
      fit2dgaus_err (columndensity, gaus2dfit, gaus2dfit_err);
      make_gaus2d_inspect (columndensity, gaus2dfit, p->shotnum.c_str ());
      gaus2d_eval_Azimuth (gaus2dfit, p->shotnum);
    }

  if (VERBOSE)
    cout << endl;

  nfit = gaus2dfit[4] * M_PI * gaus2dfit[1] * gaus2dfit[3];

  TF = p->hvbar * pow (6 * nfit, 1. / 3.);

  nfit_err =
    pow (pow (gaus2dfit_err[4] * nfit / gaus2dfit[4], 2) +
	 pow (gaus2dfit_err[1] * nfit / gaus2dfit[1],
	      2) + pow (gaus2dfit_err[3] * nfit / gaus2dfit[3], 2), 0.5);

  peakd = gaus2dfit[4] / (pow (M_PI, 0.5) * gaus2dfit[1]) / pow (p->magnif * 1e-4, 3);	// cm^-3


  if (VERBOSE)
    {
      printf ("..............  GAUSSIAN 2D FIT RESULTS ..............\n");
      printf ("ci     = %.1f pixels\n", gaus2dfit[0]);
      printf ("wi     = %.1f pixels\n", gaus2dfit[1]);
      printf ("cj     = %.1f pixels\n", gaus2dfit[2]);
      printf ("wj     = %.1f pixels\n", gaus2dfit[3]);
      printf ("peak   = %.3e \n", gaus2dfit[4]);
      printf ("offset = %.3e \n", gaus2dfit[5]);
      printf ("N from fit = %.3e \n", nfit);
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
      "------------ FIT SCATTERED PHOTONS WITH 2D GAUSISAN ------------" <<
      endl;
  if (VERBOSE)
    cout << endl;

  scatt2dfit[0] = gaus2dfit[0];
  scatt2dfit[1] = gaus2dfit[1];
  scatt2dfit[2] = gaus2dfit[2];
  scatt2dfit[3] = gaus2dfit[3];
  scatt2dfit[4] = maxNSP;
  scatt2dfit[5] = 0.1;

  if (VERBOSE)
    cout << endl << "------------ Fitting with 2D Gaussian ------------" <<
      endl;
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

  probe2dfit[0] = gaus2dfit[0];
  probe2dfit[1] = gaus2dfit[1];
  probe2dfit[2] = gaus2dfit[2];
  probe2dfit[3] = gaus2dfit[3];
  probe2dfit[4] = maxCN;

  if (VERBOSE)
    cout << endl << "------------ Fitting with 2D Gaussian ------------" <<
      endl;
  if (!p->blanks)
    fit2dgaus_no_offset (probe, probe2dfit);
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
      printf ("peak   = %.3e \n photons", probe2dfit[4]);
    }

  return;
}


void
Fermions::Fit2DFermi ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT COLUMN DENSITY WITH 2D FERMI-DIRAC ------------" <<
      endl;

  fermi2dfit[0] = gaus2dfit[4];
  fermi2dfit[1] = -5.0;
  fermi2dfit[2] = gaus2dfit[1];
  fermi2dfit[3] = gaus2dfit[3];
  fermi2dfit[4] = gaus2dfit[0];
  fermi2dfit[5] = gaus2dfit[2];
  fermi2dfit[6] = gaus2dfit[5];

  if (VERBOSE)
    {
      const char *const pnames[] =
	{ "peakd", "BetaMu", "radius_i", "radius_j", "center_i", "center_j",
	"offset"
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
      make_fermi2d_gaus2d_inspect (columndensity, fermi2dfit, gaus2dfit,
				   p->shotnum.c_str ());
      fermi2d_eval_Azimuth (fermi2dfit, p->shotnum);
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
    p->B1 / (2 * p->kbm) * ri_Fermi * ri_Fermi * p->magnif * p->magnif /
    f_Fermi;
  T_2d_ax =
    p->B2 / (2 * p->kbm) * rj_Fermi * rj_Fermi * p->magnif * p->magnif /
    f_Fermi;

  if (VERBOSE || p->showfermi)
    {
      printf ("..............  FERMI 2D FIT RESULTS ..............\n");
      printf ("n0     = %.2e\n", fermi2dfit[0]);
      printf ("BetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
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
      "------------ GET AZIMUTHAL AVERAGE OF CLOUD (ELLIPSE) ------------" <<
      endl;

  //Aspect ratio from trap geometry and self similar expansion
  double trap_aspect = sqrt (p->B2 / p->B1);

  //Aspect ratio obtained from 2D Gaussian fit
  double gaus2d_aspect = gaus2dfit[3] / gaus2dfit[1];

  // Aspect ratio determined from trap geometry and self similar expansion
  // should not be different by more than 2%
  // If they are it means that we do not know the trap frequencies
  double ar_err = abs (trap_aspect - gaus2d_aspect) / trap_aspect;
  if (ar_err > 0.02)
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
  nbins = 1024;
  binsize = 1.2;		//in pixels

  // Declare the auxiliary arrays to calculate the average 
  //    _0 is sum of values for the data
  //    _1 is number of values 
  //    average at bin =  _0 / _1
  //
  gsl_vector *azim_histo_0 = gsl_vector_alloc (nbins);
  gsl_vector *azim_histo_1 = gsl_vector_alloc (nbins);

  // Define a matrix that stores the reduced radius rho for each pixel
  gsl_matrix *rho =
    gsl_matrix_alloc (columndensity->size1, columndensity->size2);

  // Initialize arrays
  for (unsigned int j = 0; j < nbins; j++)
    {
      gsl_vector_set (azim_histo_0, j, 0.0);
      gsl_vector_set (azim_histo_1, j, 0.0);
    }

  // initialize cut arrays
  icut_r = gsl_vector_alloc (columndensity->size1);
  icut_dat = gsl_vector_alloc (columndensity->size1);

  jcut_r = gsl_vector_alloc (columndensity->size2);
  jcut_dat = gsl_vector_alloc (columndensity->size2);

  // populate auxiliary arrays for average
  // at the same time populate arrays for cuts along principal axis 
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  double xi = (double) i;	// this is y  radial
	  double xj = (double) j;	// this is g  axial

	  double data = gsl_matrix_get (columndensity, i, j);

	  // Find the radial distance, aspect ratio is taken into account
	  unsigned int dist =
	    (unsigned int)
	    floor (pow
		   (pow (p->AR * (xi - gaus2dfit[0]), 2) +
		    pow (xj - gaus2dfit[2], 2), 0.5) / binsize);


	  gsl_matrix_set (rho, i, j, dist);
	  //DEBUG: printf( " dist(%d,%d) = %d\n", i, j, dist); 

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

  string rho_path = makepath (base, p->shotnum, "_rho.TIFF");
  string rho_ascii_path = makepath (base, p->shotnum, "_rho.ascii");


  if (p->savetiff)
    toTiffImage (rho, rho_path, true);
  if (p->saveascii)
    save_gsl_matrix_ASCII (rho, rho_ascii_path);


  if (VERBOSE)
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
	  printf (" r = %.3f,  az_avg(r) = %.3f\n", r, hsum / hnum);
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


  to_dat_file_2 (azimuthal_r, azimuthal_dat, p->shotnum,
		 "datAzimuth.AZASCII");
  to_dat_file_2 (azimuthal_all_r, azimuthal_all_dat, p->shotnum,
		 "datAllAzimuth.AZASCII");
  to_dat_file_2 (icut_r, icut_dat, p->shotnum, "datIcut.AZASCII");
  to_dat_file_2 (jcut_r, jcut_dat, p->shotnum, "datJcut.AZASCII");



  return;


  /* Azimuthal-type data consists of two vector with an equan number of 
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
      "------------ FIT AZIMUTHAL AVERAGE WITH FERMI-DIRAC ------------" <<
      endl;

  /********************************* Finite temperature azimuthal fit ***/
  fermi_azimuth_fit[0] = gaus2dfit[4];
  fermi_azimuth_fit[1] = -5.0;
  //fermi_azimuth_fit[2] = pow(wi_1e*wj_1e,0.5); 
  fermi_azimuth_fit[2] = gaus2dfit[3];
  fermi_azimuth_fit[3] = gaus2dfit[5];
  fermi_azimuth_fit[4] = 0.1;


  if (VERBOSE)
    {
      const char *const pnames[] =
	{ "peakd", "BetaMu", "radius", "offset", "slope" };
      printf ("Starting fit values:\n");
      for (int e = 0; e < 5; e++)
	{
	  printf (" fermi_azimuth_fit[%d] = %f\t%s\n", e,
		  fermi_azimuth_fit[e], pnames[e]);
	}
    }

  gsl_vector *azimuthal_[2] = { azimuthal_r, azimuthal_dat };
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
      to_dat_file_2 (azimuthal_r, fitAzimuth, p->shotnum,
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
      "---------- FIT AZIMUTHAL AVERAGE WITH T=0 FERMI-DIRAC -----------" <<
      endl;
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
      gsl_vector *fitAzimuthZeroT =
	fermiAzimuthZeroT_eval (azimuthal_r, fermi_azimuth_fit_zero);
      to_dat_file_2 (azimuthal_r, fitAzimuthZeroT, p->shotnum,
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
      printf ("BetaMu     = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
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
  system (inspectstr.str ().c_str ());

  return;
}


void
Fermions::ComputeIntegrated1DDensity ()
{

  if (VERBOSE)
    cout << endl << "----------- COMPUTE INTEGRATED 1D DENSITY ------------"
      << endl;

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
      printf ("\tBetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	      fit1d_fermi_0[1], f2 (fit1d_fermi_0[1]), TF_rd);
      printf ("\tr_rd      = %.3e\n", fit1d_fermi_0[2]);
      printf ("\tc_rd      = %.3e\n", fit1d_fermi_0[3]);
      printf ("\tB_rd      = %.3e\n", fit1d_fermi_0[4]);

      printf ("\n_1 Profile:\n");
      printf ("\tn0_ax     = %.1f\n", fit1d_fermi_1[0]);
      printf ("\tBetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	      fit1d_fermi_1[1], f2 (fit1d_fermi_1[1]), TF_ax);
      printf ("\tr_ax      = %.3e\n", fit1d_fermi_1[2]);
      printf ("\tc_ax      = %.3e\n", fit1d_fermi_1[3]);
      printf ("\tB_ax      = %.3e\n", fit1d_fermi_1[4]);
    }


  return;
}
