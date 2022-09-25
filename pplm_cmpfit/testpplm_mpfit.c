/* 
 * MINPACK-1 Least Squares Fitting Library
 *
 * Test routines
 * 
 * These test routines provide examples for users to familiarize 
 * themselves with the mpfit library.  They also provide a baseline
 * test data set for users to be sure that the library is functioning
 * properly on their platform.
 *
 * By default, testmpfit is built by the distribution Makefile.
 *
 * To test the function of the mpfit library, 
 *   1. Build testmpfit   ("make testmpfit")
 *   2. Run testmpfit     ("./testmpfit")
 *   3. Compare results of your run with the distributed file testmpfit.log
 *
 * This file contains several test user functions:
 *   1. linfunc() linear fit function, y = f(x) = a - b*x
 *      - Driver is testlinfit()
 *   2. quadfunc() quadratic polynomial function, y = f(x) = a + b*x + c*x^2
 *      - Driver is testquadfit() - all parameters free
 *      - Driver is testquadfix() - linear parameter fixed
 *   3. gaussfunc() gaussian peak
 *      - Driver is testgaussfit() - all parameters free
 *      - Driver is testgaussfix() - constant & centroid fixed
 *           (this routine demonstrates in comments how to impose parameter limits)
 *   4. main() routine calls all five driver functions
 *
 * Copyright (C) 2003,2006,2009,2010, Craig Markwardt
 *
 */

/* Test routines for mpfit library
   $Id$
*/

#include "mpfit.h"
#include "pplm_mpfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

/* This is the private data structure which contains the data points
   and their uncertainties */
struct vars_struct {
  double *x;
  double *y;
  double *ey;
};

void delay(int number_of_seconds)
{
    // Converting time into milli_seconds
    int milli_seconds = 1000 * number_of_seconds;

    // Storing start time
    clock_t start_time = clock();

    // looping till required time is not achieved
    while (clock() < start_time + milli_seconds)
        ;
}
/* Simple routine to print the fit results */
void printresult(double *x, double *xact, mp_result *result) 
{
  int i;

  if ((x == 0) || (result == 0)) return;
  printf("  CHI-SQUARE = %f    (%d DOF)\n", 
	 result->bestnorm, result->nfunc-result->nfree);
  printf("        NPAR = %d\n", result->npar);
  printf("       NFREE = %d\n", result->nfree);
  printf("     NPEGGED = %d\n", result->npegged);
  printf("     NITER = %d\n", result->niter);
  printf("      NFEV = %d\n", result->nfev);
  printf("\n");
  if (xact) {
    for (i=0; i<result->npar; i++) {
      printf("  P[%d] = %f +/- %f     (ACTUAL %f)\n", 
	     i, x[i], result->xerror[i], xact[i]);
    }
  } else {
    for (i=0; i<result->npar; i++) {
      printf("  P[%d] = %f +/- %f\n", 
	     i, x[i], result->xerror[i]);
    }
  }
    
}


/* 
 * gaussian fit function
 *
 * m - number of data points
 * n - number of parameters (4)
 * p - array of fit parameters 
 *     p[0] = constant offset
 *     p[1] = peak y value
 *     p[2] = x centroid position
 *     p[3] = gaussian sigma width
 * dy - array of residuals to be returned
 * vars - private data (struct vars_struct *)
 *
 * RETURNS: error code (0 = success)
 */


int gaussfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey;
  double xc, sig2;

  x = v->x;
  y = v->y;
  ey = v->ey;

  sig2 = p[3]*p[3];

  for (i=0; i<m; i++) {
    xc = x[i]-p[2];
    dy[i] = (y[i] - p[1]*exp(-0.5*xc*xc/sig2) - p[0])/ey[i];
  }
  delay(20);
  return 0;
}



/* Test harness routine, which contains test gaussian-peak data */
int testgaussfit(int option)
{
  double x[] = {-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
		-2.8394297E-01,1.3416969E+00,1.3757038E+00,
		-1.3703436E+00,4.2581975E-02,-1.4970151E-01,
		8.2065094E-01};
  double y[] = {-4.4494256E-02,8.7324673E-01,7.4443483E-01,
		4.7631559E+00,1.7187297E-01,1.1639182E-01,
		1.5646480E+00,5.2322268E+00,4.2543168E+00,
		6.2792623E-01};
  double ey[10];
  double p[] = {0.0, 1.0, 1.0, 1.0};       /* Initial conditions */
  double pactual[] = {0.0, 4.70, 0.0, 0.5};/* Actual values used to make data*/
  double perror[4];			   /* Returned parameter errors */
  mp_par pars[4];			   /* Parameter constraints */
  int i;
  struct vars_struct v;
  int status;
  mp_result result;

  memset(&result,0,sizeof(result));      /* Zero results structure */
  result.xerror = perror;

  memset(pars,0,sizeof(pars));        /* Initialize constraint structure */
  /* No constraints */

  for (i=0; i<10; i++) ey[i] = 0.5;

  v.x = x;
  v.y = y;
  v.ey = ey;

  /* Call fitting function for 10 data points and 4 parameters (no
     parameters fixed) */

  if (option==0) {
     printf("CALLING mpfit()\n");
     status = mpfit(gaussfunc, 10, 4, p, pars, 0, (void *) &v, &result);
     printf("*** testgaussfit (mpfit) status = %d\n", status);
  } else {
     printf("CALLING pplm_mpfit()\n");
     status = pplm_mpfit(gaussfunc, 10, 4, p, pars, 0, (void *) &v, &result);
     printf("*** testgaussfit (pplm_mpfit) status = %d\n", status);
  }
  printresult(p, pactual, &result);

  return 0;
}




/* Test harness routine, which contains test gaussian-peak data */
int test_pplm()
{
  double x[] = {-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
		-2.8394297E-01,1.3416969E+00,1.3757038E+00,
		-1.3703436E+00,4.2581975E-02,-1.4970151E-01,
		8.2065094E-01};
  double y[] = {-4.4494256E-02,8.7324673E-01,7.4443483E-01,
		4.7631559E+00,1.7187297E-01,1.1639182E-01,
		1.5646480E+00,5.2322268E+00,4.2543168E+00,
		6.2792623E-01};
  double ey[10];
  double p[] = {0.0, 1.0, 1.0, 1.0};       /* Initial conditions */
  double pactual[] = {0.0, 4.70, 0.0, 0.5};/* Actual values used to make data*/
  double perror[4];			   /* Returned parameter errors */
  mp_par pars[4];			   /* Parameter constraints */
  int i;
  struct vars_struct v;
  int status;
  mp_result result;

  memset(&result,0,sizeof(result));      /* Zero results structure */
  result.xerror = perror;

  memset(pars,0,sizeof(pars));        /* Initialize constraint structure */
  /* No constraints */

  for (i=0; i<10; i++) ey[i] = 0.5;

  v.x = x;
  v.y = y;
  v.ey = ey;

  /* Call fitting function for 10 data points and 4 parameters (no
     parameters fixed) */
  status = pplm_mpfit(gaussfunc, 10, 4, p, pars, 0, (void *) &v, &result);

  printf("*** test PPWLM gaussfit status = %d\n", status);
  printresult(p, pactual, &result);

  return 0;
}


/* Main function which drives the whole thing */
int main(int argc, char *argv[])
{
  float time0, time1,time2;
  time0=(float)clock()/CLOCKS_PER_SEC;

  testgaussfit(0);
  time1=(float)clock()/CLOCKS_PER_SEC;
  printf("Run Time=%f\n",time1-time0);

  testgaussfit(1);
  time2=(float)clock()/CLOCKS_PER_SEC;
  printf("Run Time=%f\n",time2-time1);

  exit(0);
}
