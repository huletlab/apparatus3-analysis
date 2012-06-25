/*
 *
 *  Evaluator class for fluorescence pictuers taken with the Basler camera in Apparatus 3
 *
 *
 *  Author : Pedro M Duarte  Jul 2011
 *
 * 
 */

// This structure contains all parameters relevant for analysis

// Convention is  size1 : 0 : radial : i
//                size2 : 1 : axial  : j 

extern bool VERBOSE;


struct params
{
  unsigned int shot;

  string shotnum, reportfile, atomsfile;

  unsigned int roi[4], roisize[2], centerpt[2];

  bool verbose, center, crop, plots, reanalyze, roi_user, roisize_user;

  double lambda, hc, isat, gamma, tau, quantumeff, solidangle, magnif, texp, det, detcalib;

};



class Basler
{
public:

  Basler (struct params *params)
  {
    p = params;
  }
   ~Basler ()
  {
    return;
  }

  void LoadFLUOR ();
  void ComputeColumnDensity ();
  void NAtoms ();
  double GetNPixels ()
  {
    return atoms->size1 * atoms->size2;
  }
  void Fit1DGauss (bool col_row);
  void Fit2DGauss ();
  void ComputeRadialAxialDensity ();

  double number, number_fit, maxOD, nsp, gaus2dfit[6], gaus2dfit_err[6];
  double abs_ci, abs_cj;	// centers of cloud in the uncropped pict

  struct params *p;

private:
  gsl_matrix * atoms;
  gsl_matrix * columndensity;
  gsl_matrix *residuals;

  gsl_vector *axialdensity;
  gsl_vector *radialdensity;

  unsigned int ci, cj, FWHMi, FWHMj;
  double offset, peak, ci_, cj_, wi_1e, wj_1e;
};


void
Basler::LoadFLUOR ()
{
  abs_ci = 0.;
  abs_cj = 0.;

  atoms = ReadFluorImg_gsl_matrix (p->atomsfile);

  if (p->roisize_user)
    {

      if (p->roisize[0] > atoms->size1 || p->roisize[1] > atoms->size2)
	{
	  cout << "\tERROR:  Size of ROI is larger than image" << endl;
	  exit (2);
	}

      p->roi[2] = p->roisize[0];
      p->roi[3] = p->roisize[1];
      ComputeColumnDensity ();
      if (VERBOSE)
	{
	  cout << endl << "------------ Finding ROI box ------------";
	}
      findpeak (columndensity, &ci, &cj, &peak);
      gsl_matrix_free (columndensity);

      if (VERBOSE)
	{
	  cout << "\tROI Size = " << p->roi[2] << ", " << p->roi[3] << endl;
	  cout << "\tci=" << ci << ",  cj=" << cj << endl;
	  cout << "\tSTART BOX = [ " << (int) ci -
	    (int) p->roisize[0] / 2 << " : " << (int) ci +
	    (int) p->roisize[0] / 2 << " , " << (int) cj -
	    (int) p->roisize[1] / 2 << " : " << (int) cj +
	    (int) p->roisize[1] / 2 << " ]" << endl;
	}

      while ((int) ci - (int) p->roi[2] / 2 > 0
	     && (int) ci + (int) p->roi[2] / 2 > (int) atoms->size1)
	{
	  ci--;
	}

      while ((int) cj - (int) p->roi[3] / 2 > 0
	     && (int) cj + (int) p->roi[3] / 2 > (int) atoms->size2)
	{
	  cj--;
	}

      while ((int) ci - (int) p->roi[2] / 2 < 0
	     && (int) ci + (int) p->roi[2] / 2 < (int) atoms->size1)
	{
	  ci++;
	}

      while ((int) cj - (int) p->roi[3] / 2 < 0
	     && (int) cj + (int) p->roi[3] / 2 < (int) atoms->size2)
	{
	  cj++;
	}

      if (VERBOSE || true)
	{
	  cout << "\tFINAL BOX = [ " << (int) ci -
	    (int) p->roisize[0] / 2 << " : " << (int) ci +
	    (int) p->roisize[0] / 2 << " , " << (int) cj -
	    (int) p->roisize[1] / 2 << " : " << (int) cj +
	    (int) p->roisize[1] / 2 << " ]" << endl;
	}

      p->roi[0] = ci - p->roisize[0] / 2;
      p->roi[1] = cj - p->roisize[1] / 2;

      p->roi_user = true;
      if (VERBOSE)
	{
	  cout << "\tROI found at " << p->roi[0] << "," << p->
	    roi[1] << "," << p->roi[2] << "," << p->roi[3] << endl;
	}

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
      gsl_matrix_free (atoms);
      atoms = catoms;
    }


  if (VERBOSE)
    cout << endl;

  return;
}


void
Basler::ComputeColumnDensity ()
{

  unsigned int s1 = atoms->size1;
  unsigned int s2 = atoms->size2;


  columndensity = gsl_matrix_alloc (s1, s2);
  
  double pee=0.5; 

  double c2a =  p->detcalib*(p->tau/ p->texp)*(p->quantumeff/p->solidangle) / pee; 
  
  double at; 

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{

	  at = gsl_matrix_get (atoms, i, j);  //counts 
	  gsl_matrix_set (columndensity, i, j, at*c2a);  //convert counts -> atoms

	}
    }

  if (VERBOSE)
    {
      cout << endl << "------------ Column Density Stats ------------" <<
	endl;
      cout << "\ts1=" << s1 << ", s2=" << s2 << endl;
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

  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string column_path ("");
  column_path += path;
  column_path += "/";
  column_path += "fluor_";
  column_path += p->shotnum;
  string column_ascii_path ("");
  column_ascii_path = column_path;
  column_path += ".TIFF";
  column_ascii_path += ".ascii";
  toTiffImage (columndensity, column_path);
  save_gsl_matrix_ASCII (columndensity, column_ascii_path);
  return;
}


void
Basler::Fit1DGauss (bool col_row)
{


  //True fits a column , False fits a row
  unsigned int ii = col_row ? cj : ci;
  if (VERBOSE)
    {
      cout << endl << "------------ Fitting 1D Gaussian ------------" << endl;
      cout << "Column density ( " << columndensity->
	size1 << ", " << columndensity->size2 << ")  ";
      if (col_row)
	cout << "Fitting COLUMN " << ii << endl;
      if (!col_row)
	cout << "Fitting ROW " << ii << endl;
    }

  if (col_row && ii > columndensity->size2 || ii < 0)
    {
      cout << "Cannot fit column with index " << ii << " : out of bounds!" <<
	endl;
      return;
    }
  if (!col_row && ii > columndensity->size1 || ii < 0)
    {
      cout << "Cannot fit row with index " << ii << " : out of bounds!" <<
	endl;
      return;
    }

  unsigned int s = col_row ? columndensity->size1 : columndensity->size2;

  gsl_vector *profile = gsl_vector_alloc (s);
  for (unsigned int index = 0; index < s; index++)
    {
      unsigned int i = col_row ? index : ii;
      unsigned int j = col_row ? ii : index;
      gsl_vector_set (profile, index, gsl_matrix_get (columndensity, i, j));
    }

  //findpeak (columndensity, &ci, &cj, &peak);
  //findFWHM (columndensity, &FWHMi, &FWHMj);

  double center = (double) (col_row ? ci : cj);
  double FWHM = (double) (col_row ? FWHMi : FWHMj);


  if (VERBOSE)
    cout << endl;
  double gausfit1d[4] = { center, FWHM, peak, 0.1 };
  fit1dgaus (profile, gausfit1d);

  gausfit1d[1] = fabs (gausfit1d[1]);

  if (col_row)
    {
      ci_ = gausfit1d[0];
      wi_1e = gausfit1d[1];
    }
  if (!col_row)
    {
      cj_ = gausfit1d[0];
      wj_1e = gausfit1d[1];
    }

  gsl_vector_free (profile);

  return;
}



void
Basler::Fit2DGauss ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT COLUMN DENSITY WITH 2D GAUSISAN ------------" << endl;
  //findpeak (columndensity, &ci, &cj, &peak);
  findcenter (columndensity, &ci, &cj, &peak);
  findFWHM (columndensity, &FWHMi, &FWHMj);

  if (VERBOSE)
    cout << endl;
  Fit1DGauss (false);		//False is to fit row
  if (VERBOSE)
    cout << endl;
  Fit1DGauss (true);		//True is to fit column

  gaus2dfit[0] = ci_;
  gaus2dfit[1] = wi_1e;
  gaus2dfit[2] = cj_;
  gaus2dfit[3] = wj_1e;
  gaus2dfit[4] = peak;
  gaus2dfit[5] = 0.1;



  if (VERBOSE)
    cout << endl << "------------ Fitting with 2D Gaussian ------------" <<
      endl;
  fit2dgaus_err (columndensity, gaus2dfit, gaus2dfit_err);
  if (VERBOSE)
    cout << endl;

  gaus2dfit[1] = fabs (gaus2dfit[1]);
  gaus2dfit[3] = fabs (gaus2dfit[3]);

  ci_ = gaus2dfit[0];
  wi_1e = gaus2dfit[1];
  cj_ = gaus2dfit[2];
  wj_1e = gaus2dfit[3];
  peak = gaus2dfit[4];
  offset = gaus2dfit[5];

  abs_ci += ci_;
  abs_cj += cj_;

  residuals = gsl_matrix_alloc (columndensity->size1, columndensity->size2);

  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);


  FILE *fitI, *fitJ;
  string fitI_path (path);
  fitI_path += "/";
  fitI_path += p->shotnum;
  fitI_path += "_fitI.fitdat";
  string fitJ_path (path);
  fitJ_path += "/";
  fitJ_path += p->shotnum;
  fitJ_path += "_fitJ.fitdat";

  fitI = fopen (fitI_path.c_str (), "w+");
  fitJ = fopen (fitJ_path.c_str (), "w+");

  fprintf (fitI, "#pixel\tdat\tmodel\n");
  fprintf (fitJ, "#pixel\tdat\tmodel\n");

  ci = (unsigned int) floor (ci_) - 1;
  cj = (unsigned int) floor (cj_) - 1;

  ci = coerce_matrix_index (ci, columndensity->size1);
  cj = coerce_matrix_index (cj, columndensity->size2);

  if (VERBOSE)
    {
      cout << "ci = " << ci << " ; cj = " << cj << endl;
    }

  number_fit = 0.0;

  double model, data;
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  double xi = (double) i;
	  double xj = (double) j;
	  model =
	    offset +
	    peak * exp (-1. *
			(pow ((xi - ci_) / wi_1e, 2.) +
			 pow ((xj - cj_) / wj_1e, 2.)));
	  data = gsl_matrix_get (columndensity, i, j);
	  gsl_matrix_set (residuals, i, j, data - model);
	  if (i == ci)
	    fprintf (fitJ, "%e\t%e\t%e\n", xj, data, model);
	  if (j == cj)
	    fprintf (fitI, "%e\t%e\t%e\n", xi, data, model);

	  number_fit += model;

	}
    }

  if (VERBOSE)
    {
      cout << endl << "number_fit = " << number_fit << endl;
    }

  fclose (fitI);
  fclose (fitJ);

  std::ofstream gpl;
  gpl.open ("temp.gpl");
  gpl << "set size 1.0,0.45" << endl;
  gpl << "set multiplot" << endl;
  gpl << "set origin 0.0,0.0" << endl;
  gpl << "plot \"" << p->
    shotnum << "_fitI.fitdat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_fitI.fitdat\" u 1:3 w lines title \"" << p->
    shotnum << "FIT I\" " << endl;
  gpl << "set origin 0.0,0.5" << endl;
  gpl << "plot \"" << p->
    shotnum << "_fitJ.fitdat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_fitJ.fitdat\" u 1:3 w lines title \"" << p->
    shotnum << "FIT J\" " << endl;
  gpl << "unset multiplot" << endl;
  if (p->plots)
    std::system ("gnuplot -persist temp.gpl");
  gpl.close ();
  remove ("temp.gpl");

  gpl.open ("temp.gpl");
  gpl << "set terminal png medium" << endl;
  gpl << "set output \"cuts_" << p->shotnum << ".png\"" << endl;
  gpl << "set size 1.0,1.0" << endl;
  gpl << "set multiplot" << endl;
  gpl << "set origin 0.0,0.0" << endl;
  gpl << "set size 1.0,0.45" << endl;
  gpl << "plot \"" << p->
    shotnum << "_fitI.fitdat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_fitI.fitdat\" u 1:3 w lines title \"" << p->
    shotnum << "FIT I\" " << endl;
  gpl << "set origin 0.0,0.5" << endl;
  gpl << "set size 1.0,0.45" << endl;
  gpl << "plot \"" << p->
    shotnum << "_fitJ.fitdat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_fitJ.fitdat\" u 1:3 w lines title \"" << p->
    shotnum << "FIT J\" " << endl;
  gpl << "unset multiplot" << endl;
  gpl.close ();
  std::system ("gnuplot temp.gpl");
  remove ("temp.gpl");

  string residuals_path ("");
  residuals_path += path;
  residuals_path += "/";
  residuals_path += "res_";
  residuals_path += p->shotnum;
  residuals_path += ".TIFF";
  toTiffImage (residuals, residuals_path);


  return;
}

void
Basler::ComputeRadialAxialDensity ()
{

  unsigned int s1 = columndensity->size1;
  unsigned int s2 = columndensity->size2;

  radialdensity = gsl_vector_alloc (s1);
  axialdensity = gsl_vector_alloc (s2);

  gsl_vector_set_all (radialdensity, 0.0);
  gsl_vector_set_all (axialdensity, 0.0);

  double cd_ij = 0.0;

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  cd_ij = gsl_matrix_get (columndensity, i, j);

	  gsl_vector_set (radialdensity, i,
			  gsl_vector_get (radialdensity, i) + cd_ij);
	  gsl_vector_set (axialdensity, j,
			  gsl_vector_get (axialdensity, j) + cd_ij);
	}
    }

  if (VERBOSE)
    cout << endl;

  double radialfit[4] = { ci, wi_1e, gsl_vector_max (radialdensity), 0.1 };
  fit1dgaus (radialdensity, radialfit);

  double axialfit[4] = { cj, wj_1e, gsl_vector_max (axialdensity), 0.1 };
  fit1dgaus (axialdensity, axialfit);


  if (VERBOSE)
    {
      cout << endl <<
	"------------ Radial and Axial Density Stats ------------" << endl;
      cout << "\tc_rd=" << radialfit[0] << ", w_rd=" << radialfit[1] << ", A="
	<< radialfit[2] << ", B=" << radialfit[3] << endl;
      cout << "\tc_ax=" << axialfit[0] << ", w_ax=" << axialfit[1] << ", A="
	<< axialfit[2] << ", B=" << axialfit[3] << endl;
    }



  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);

  FILE *radial_F, *axial_F;
  string radial_path (path);
  radial_path += "/";
  radial_path += p->shotnum;
  radial_path += "_radial.ndat";
  string axial_path (path);
  axial_path += "/";
  axial_path += p->shotnum;
  axial_path += "_axial.ndat";

  radial_F = fopen (radial_path.c_str (), "w+");
  axial_F = fopen (axial_path.c_str (), "w+");

  fprintf (radial_F, "#pixel\tdat\tmodel\n");
  fprintf (axial_F, "#pixel\tdat\tmodel\n");

  double model, data;

  for (unsigned int i = 0; i < s1; i++)
    {
      double xi = (double) i;
      model =
	radialfit[3] +
	radialfit[2] * exp (-1.0 *
			    pow ((xi - radialfit[0]) / radialfit[1], 2.));
      data = gsl_vector_get (radialdensity, i);
      fprintf (radial_F, "%e\t%e\t%e\n", xi, data, model);
    }

  for (unsigned int j = 0; j < s2; j++)
    {
      double xj = (double) j;
      model =
	axialfit[3] +
	axialfit[2] * exp (-1.0 * pow ((xj - axialfit[0]) / axialfit[1], 2.));
      data = gsl_vector_get (axialdensity, j);
      fprintf (axial_F, "%e\t%e\t%e\n", xj, data, model);
    }

  fclose (radial_F);
  fclose (axial_F);

  // Save plots to png file
  std::ofstream gpl;
  gpl.open ("temp.gpl");
  gpl << "set terminal png medium" << endl;
  gpl << "set output \"radial_" << p->shotnum << ".png\"" << endl;
  gpl << "set size 1.0,0.45" << endl;
  gpl << "plot \"" << p->
    shotnum << "_radial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_radial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << " radial\" " << endl;
  gpl << "set output \"axial_" << p->shotnum << ".png\"" << endl;
  gpl << "plot \"" << p->
    shotnum << "_axial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_axial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << " axial\" " << endl;
  std::system ("gnuplot -persist temp.gpl");
  gpl.close ();

  // Show plots on screen if the plots option is selected 
  gpl.open ("temp.gpl");
  gpl << "set size 1.0,0.45" << endl;
  gpl << "set multiplot" << endl;
  gpl << "set origin 0.0,0.0" << endl;
  gpl << "plot \"" << p->
    shotnum << "_radial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_radial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << "radial\" " << endl;
  gpl << "set origin 0.0,0.5" << endl;
  gpl << "plot \"" << p->
    shotnum << "_axial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_axial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << "axial\" " << endl;
  gpl << "unset multiplot" << endl;
  gpl.close ();
  if (p->plots)
    std::system ("gnuplot -persist temp.gpl");

  return;
}
