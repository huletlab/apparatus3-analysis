
#include <omp.h>
#include <math.h>

double
fermi2d_simplex_f (const gsl_vector * v, void *params)
{

  //benchmark start
//  double start = omp_get_wtime(); 
  
  unsigned int s1 = ((gsl_matrix *) params)->size1;
  unsigned int s2 = ((gsl_matrix *) params)->size2;

  double n0     = gsl_vector_get (v, 0);
  double BetaMu = gsl_vector_get (v, 1);
  double ri     = gsl_vector_get (v, 2);
  double rj     = gsl_vector_get (v, 3);
  double ci     = gsl_vector_get (v, 4);
  double cj     = gsl_vector_get (v, 5);
  double B      = gsl_vector_get (v, 6); 

  double sumsq = 0.;

  //i is radial === y
  //j is axial  === g

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  double model = n0/ f1(BetaMu) * f1( BetaMu - fq(BetaMu) * ( pow ( (j -cj)/rj ,2) + pow( ( i-ci)/ri ,2))) + B; 
	  double dat = gsl_matrix_get ((gsl_matrix *) params, i, j);
	  sumsq += pow (model - dat, 2);
	}
    }

  //end benchmark
//  double end = omp_get_wtime();
//  printf( "time elapsed in fermi errfunc %.5f\n", end-start); 

  return sumsq;
}

double
fermi1d_azimuthal_simplex_f (const gsl_vector * v,  void *params)
{


  //benchmark start
//  double start = omp_get_wtime(); 
  
  unsigned int s = (((gsl_vector **) params)[0])->size; 
  

  double n0     = gsl_vector_get (v, 0);
  double BetaMu = gsl_vector_get (v, 1);
  double r     = gsl_vector_get (v, 2);
  double B      = gsl_vector_get (v, 3);
  double mx    = gsl_vector_get(v, 4); 

  double sumsq = 0.;
  
  gsl_vector *d  = ((gsl_vector **)params)[0]; 
  gsl_vector *az = ((gsl_vector **)params)[1]; 

  //i is radial === y
  //j is axial  === g

  for (unsigned int i = 0; i < s; i++){
  //printf("Inisde simplex_f loop\n"); 
          double dist = gsl_vector_get(d,i); 
	  double dat = gsl_vector_get (az, i);
	  double model = n0/ f1(BetaMu) * f1( BetaMu - fq(BetaMu) * ( pow ( dist/r ,2) )) + B + mx*dist; 
	  sumsq += pow (model - dat, 2);
    }

//  end benchmark
//  double end = omp_get_wtime();
//  printf( "time elapsed in fermi errfunc %.5f\n", end-start); 

  return sumsq;
}

double
fermi1d_azimuthal_zero_simplex_f (const gsl_vector * v,  void *params)
{
  //benchmark start
//  double start = omp_get_wtime(); 
  
  unsigned int s = (((gsl_vector **) params)[0])->size; 
  

  double n0     = gsl_vector_get (v, 0);
  double r     = gsl_vector_get (v, 1);
  double B      = gsl_vector_get (v, 2);
  double mx    = gsl_vector_get(v, 3); 

  double sumsq = 0.;
  
  gsl_vector *d  = ((gsl_vector **)params)[0]; 
  gsl_vector *az = ((gsl_vector **)params)[1]; 

  //i is radial === y
  //j is axial  === g

  for (unsigned int i = 0; i < s; i++){
  //printf("Inisde simplex_f loop\n"); 
          double dist = gsl_vector_get(d,i); 
	  double dat = gsl_vector_get (az, i);
          double model = n0 * pow( std::max ( 1. - pow( dist/r,2.) , 0.), 2.) + B + mx*dist;    
	  sumsq += pow (model - dat, 2);
    }

//  end benchmark
//  double end = omp_get_wtime();
//  printf( "time elapsed in fermi errfunc %.5f\n", end-start); 

  return sumsq;
}

double
fermi1d_simplex_f (const gsl_vector *v, void *params)
{
  unsigned int s1 = ((gsl_vector *) params)->size;

  double n0     = gsl_vector_get (v, 0);
  double BetaMu = gsl_vector_get (v, 1);
  double r      = gsl_vector_get (v, 2);
  double c      = gsl_vector_get (v, 3);
  double B      = gsl_vector_get (v, 4); 

  double sumsq = 0.;

  //i is radial === y
  //j is axial  === g

  for (unsigned int i = 0; i < s1; i++)
    {
	  double model = n0/ f32(BetaMu) * f32( BetaMu - fq(BetaMu) *  pow ( (i -c)/r ,2) ) + B ; 
	  double dat = gsl_vector_get ((gsl_vector *) params, i);
	  sumsq += pow (model - dat, 2);
    }

  return sumsq;
}



/*

int
fermi2d_f (const gsl_vector * x, void *data, gsl_vector * f)
{

  // x contains the parameters
  // data  contains the column density matrix
  // f is a vector containing the differences between matrix and fit function

  unsigned int s1 = ((gsl_matrix *) data)->size1;
  unsigned int s2 = ((gsl_matrix *) data)->size2;

  // Trap parameters 
  wx = 2.*M_PI*3800.;
  wy = 2.*M_PI*3800.;
  wz = 2.*M_PI*3800./8.;
  a  = cos(52.5*M_PI/180.);
  b  = sin(52.5*M_PI/180.);

  bx = pow(wx,2)/(1+pow(wx*t,2));
  by = pow(wy,2)/(1+pow(wy*t,2));
  bz = pow(wz,2)/(1+pow(wz*t,2));

  magnif = 3.2;  // um per pixel

  A  = m*( pow(by,0.5) *( pow(a,2)*bx + pow(b,2)*bz - pow(a*b,2)*pow(bz-bx,2)/(pow(b,2)*bx + pow(a,2)*bz)) );
  B1 = m*( pow(a,2)*bx + pow(b,2)*bz - pow(a*b,2)*pow(bz-bx,2)/(pow(b,2)*bx + pow(a,2)*bz));
  B2 = m*( bz );

  double N      = gsl_vector_get(x,0);
  double Beta   = gsl_vector_get(x,1);
  double BetaMu = gsl_vector_get(x,2);
  double cy     = gsl_vector_get(x,3);
  double cg     = gsl_vector_get(x,4);

  size_t ii = 0;

  //i is radial === y
  //j is axial  === g

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
          double model = 
            N*Beta/(2.*M_PI*F2(BetaMu)) * A * F1(BetaMu - Beta/2.* (  B1*pow( (j-cg)*magnif , 2) + B2*pow( (i-cy)*magnif , 2))); 
	  double dat = gsl_matrix_get ((gsl_matrix *) data, i, j);
	  gsl_vector_set (f, ii, (model - dat));
	  ii++;
	}
    }
}


*/
