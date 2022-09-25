/*
 ====================================================
   Parallel Processing of Levenberg Marquard Methods
             via Wrapping functions
 ====================================================
 */
#include "mpfit.h"
#include <math.h>
#include <stdlib.h>

/* Variable Local To these PPLM functions */
mp_func *pplm_kernel_func_ptr;/*Adress of the given fitting function */
double  *pplm_p_last;         /*Adress of values of the last par vals*/
double  *pplm_d_last;         /*Adress of values of the last return vals*/
int      pplm_ncalls;         /*How many times pplm_wrap has been called*/

int pplm_wrap(mp_func funct, int m, int n, double *p, double *dy,
                                           double **dvec, void *vars) {

  static double zero = 0.0;
  static double **nptr=0;
  static double eps = sqrt(MP_MACHEP0);
  static int    f_last;
  int iflag[n];
  double sum=0.0;

  /* Increment the call count of this subroutine since pplm_mpfit
   * initialisation*/
  pplm_ncalls++;

  /* If not the first call check if given parameter value set is the
   *  same as those in the last call */
  if (pplm_ncalls>1) {
    sum=0.0;
    for (int i=0;i<n;i++) sum=sum+pow(p[i]-pplm_p_last[i],2);
  }

  if (pplm_ncalls>1 && sum==0.0) {
    /* The input parameters are the same as those used in the last call
       so we copy and return the last result values */
    iflag[0]=f_last;
    for (int i=0;i<m;i++) dy[i]=pplm_d_last[i];
  } else {
    /* The input parameters are not the same as those used in the last call
       so we call the object function to derive new result values */
    iflag[0]=funct(m,n,p,dy,dvec,vars);

    /* Store current parameter values and the results for the next call */
    for (int i=0;i<n;i++) pplm_p_last[i]=p[i];
    for (int i=0;i<m;i++) pplm_d_last[i]=dy[i];
    f_last=iflag[0];
  }

  /* If dvec is null then Jacobian elements are not required so return*/
  if (!dvec ) return iflag[0];

  #pragma omp parallel for
  for (int j=0;j<n;j++) {
    double pl[n];
    for (int i=0; i<n; i++) pl[i]=p[i];
    double temp = pl[j];
    double h = eps * fabs(temp);
    if (h == zero) h = eps;
    pl[j]=pl[j]+h;
    iflag[j]=funct(m,n,pl,&dvec[j][0],nptr,vars);
    for (int i=0; i<m;i++) dvec[j][i]=(dvec[j][i]-dy[i])/h;
  }

  /* In the event of failure return the first failed flag*/
  for (int j=0; j<n; j++) if (iflag[j]!=0) return iflag[j];

  return 0;
}

int pplm_shell_func(int m, int n, double *p, double *dy,
                                           double **dvec, void *vars) {

  pplm_wrap(*pplm_kernel_func_ptr,m,n,p,dy,dvec,vars);

}

int pplm_mpfit(mp_func funct, int m, int npar,
	  double *xall, mp_par *pars, mp_config *config, void *private_data,
	  mp_result *result) {

  int iflag;

  /* Assign the private variables static to pplm */
  pplm_kernel_func_ptr = &funct;
  pplm_p_last=(double *)malloc(npar*sizeof(double));
  pplm_d_last=(double *)malloc(   m*sizeof(double));
  pplm_ncalls=0;

  mp_par new_par_config[npar];

  /* Make a copy of the par configurations so that side==3
   * (i.e. tells mpfit that all Jacobian variables are analytically
   * derived) */
  for (int i=0;i<npar;i++) {
    new_par_config[i].fixed        = pars[i].fixed;
    new_par_config[i].limited[0]   = pars[i].limited[0];
    new_par_config[i].limited[1]   = pars[i].limited[1];
    new_par_config[i].limits[0]    = pars[i].limits[0];
    new_par_config[i].limits[1]    = pars[i].limits[1];
    new_par_config[i].parname      = "/0";
    new_par_config[i].step         = pars[i].step;
    new_par_config[i].side         = 3;
    new_par_config[i].deriv_debug  = 0;
    new_par_config[i].deriv_reltol = pars[i].deriv_reltol;
    new_par_config[i].deriv_abstol = pars[i].deriv_abstol;
  }

  /* Call mpfit on the pplm_shell_func */
  iflag=mpfit(pplm_shell_func,m,npar,xall,new_par_config,config,private_data,result);

  /* Cleanup memory allocations */
  free(pplm_p_last);
  free(pplm_d_last);

  /* Return the mpfit output value */
  return iflag;

}


