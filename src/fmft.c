/* This program implements the Frequency Modified Fourier Transform 
(Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137). 
Given a quasi--periodic complex signal X + iY, the algorithm 
estimates the frequencies (f_j), amplitudes (A_j) and phases 
(psi_j) in its decomposition:

X(t) + iY(t) = Sum_j=1^N [ A_j * exp i (f_j * t + psi_j) ] */      

#define FMFT_TOL 1.0e-10 /* MFT NOMINAL PRECISION */
#define FMFT_NEAR 0.     /* MFT OVERLAP EXCLUSION PARAMETER */

#include <stdio.h>
#include <math.h>

#include "nrutil.h"
#include "stdlib.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h> 

#define PI 3.14159265358979
#define TWOPI (2.*PI)


static int itemp;
static unsigned long ultemp;
static float ftemp;
static double dtemp;

#define FSQR(a) ((ftemp=(a)) == 0.0 ? 0.0 : ftemp*ftemp)
#define DSQR(a) ((dtemp=(a)) == 0.0 ? 0.0 : dtemp*dtemp)

#define SHFT3(a,b,c) (a)=(b);(b)=(c)
#define SHFT4(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

#define ISWAP(a,b) itemp=(a);(a)=(b);(b)=itemp
#define ULSWAP(a,b) ultemp=(a);(a)=(b);(b)=ultemp
#define FSWAP(a,b) ftemp=(a);(a)=(b);(b)=ftemp

void window(double *x, double *y, double *xdata, double *ydata, size_t ndata);

void power(float *powsd, double *x, double *y, size_t ndata);

double bracket(float *powsd, size_t ndata);

double golden(double (*f)(double, double *, double *, size_t), 
	      double leftf, double centerf, double rightf, 
	      double *x, double *y, size_t ndata);

void phifun(double *xphi, double *yphi, double freq,  
	      double xdata[], double ydata[], long n);

double phisqr(double freq, double xdata[], double ydata[], size_t ndata);

void amph(double *amp, double *phase, double freq, 
	  double xdata[], double ydata[], size_t ndata);


void dindex(unsigned long n, double arr[], unsigned long indx[]);



/* THE MAIN FUNCTION ****************************************************/

  struct component{
    double freq;
    double amp;
    double phase;
  };
 

int fmft(int *localnfreq, double *localminfreq, double *localmaxfreq, int *localflag, 
	 int *localndata, double *localxdata, double *localydata,  
   struct component *signal1, struct component *signal2, struct component *signal3);


int fmft(int *localnfreq, double *localminfreq, double *localmaxfreq, int *localflag, 
	 int *localndata, double *localxdata, double *localydata,  
   struct component *signal1, struct component *signal2, struct component *signal3)

/*  struct component signal1[nfreq];*/
/*  struct component signal2[nfreq];*/
/*  struct component signal3[nfreq];*/
/*   struct component* signal1  = new struct component[localnfreq];*/
/*   struct component* signal2  = new struct component[localnfreq];*/
/*   struct component* signal3  = new struct component[localnfreq];*/

/* 
MC : signal1, signal2 and signal3 are now replacing the old output array
In the output array **output: output[3*flag-2][i], output[3*flag-1][i] 
and output[3*flag][i] are the i-th frequency, amplitude and phase; nfreq is the 
number of frequencies to be computed (the units are rad/sep, where sep is the 
`time' separation between i and i+1. The algorithm is  

Basic Fourier Transform algorithm           if   flag = 0;   not implemented   
Modified Fourier Transform                  if   flag = 1;
Frequency Modified Fourier Transform        if   flag = 2;
FMFT with additional non-linear correction  if   flag = 3

(while the first algorithm is app. 3 times faster than the third one, 
the third algorithm should be in general much more precise).  
The computed frequencies are in the range given by minfreq and maxfreq.
The function returns the number of determined frequencies or 0 in the case
of error.

The vectors input[1][j] and input[2][j], j = 1 ... ndata (ndata must
be a power of 2), are the input data X(j-1) and Y(j-1).
*/   
     
{
  int nearfreqflag;
  long i,j,k,l,m;
  float *powsd;
  double *xdata, *ydata, *x, *y;
  double centerf, leftf, rightf, fac, xsum, ysum;
  double **freq, **amp, **phase, *f, *A, *psi;
  double **Q, **alpha, *B;

  FILE *fp;

  int nfreq = *localnfreq;
  double minfreq = *localminfreq;
  double maxfreq = *localmaxfreq;
  int flag = *localflag;
  size_t ndata = *localndata;

/*  printf("prelimarg %d, %zu %d", nfreq, ndata, flag);*/
  printf("prelimarg %d, %d %d", nfreq, *localndata, *localflag);
  printf("xdata: %.2f, %.2f, %.2f", localxdata[0],localxdata[1],localxdata[2]);
  printf("ydata: %.2f, %.2f, %.2f", localydata[0],localydata[1],localydata[2]);

  /* ALLOCATION OF VARIABLES */

/*  xdata = dvector(1,ndata);*/
/*  ydata = dvector(1,ndata);*/
  x = dvector(1,ndata);
  y = dvector(1,ndata);

  powsd = vector(1, ndata);
  
  freq = dmatrix(1, 3*flag, 1, nfreq); 
  amp = dmatrix(1, 3*flag, 1, nfreq);
  phase = dmatrix(1, 3*flag, 1, nfreq);

  f = dvector(1, nfreq);
  A = dvector(1, nfreq);
  psi = dvector(1, nfreq);

  
  Q = dmatrix(1, nfreq, 1, nfreq); 
  alpha = dmatrix(1, nfreq, 1, nfreq);
  B = dvector(1, nfreq);
 

  /* 1 LOOP FOR MFT, 2 LOOPS FOR FMFT, 3 LOOPS FOR NON-LINEAR FMFT */

  for(l=1; l<=flag; l++){
 
    if(l==1){
      xdata = localxdata -1;  // -1 because dvector vs *double
      ydata = localydata -1;
      /* SEPARATE REAL AND IMAGINERY PARTS */ 
/*      for(j=1;j<=ndata;j++){*/
/*  xdata[j] = localxdata[j-1];*/
/*  ydata[j] = localydata[j-1];*/

    } else {

       /* GENERATE THE QUASIPERIODIC FUNCTION COMPUTED BY MFT */
      for(i=1;i<=ndata;i++){
	xdata[i] = 0; ydata[i] = 0; 
	for(k=1;k<=nfreq;k++){
	  xdata[i] += amp[l-1][k]*cos(freq[l-1][k]*(i-1) + phase[l-1][k]);
	  ydata[i] += amp[l-1][k]*sin(freq[l-1][k]*(i-1) + phase[l-1][k]);
	}
      }

    }
  

    /* MULTIPLY THE SIGNAL BY A WINDOW FUNCTION, STORE RESULT IN x AND y */
    window(x, y, xdata, ydata, ndata);
    
    /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
    power(powsd, x, y, ndata);


    if(l==1) 

      /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
      while((centerf = bracket(powsd, ndata)) < minfreq || centerf > maxfreq) {

	
	/* IF NO, SUBSTRACT IT FROM THE SIGNAL */
	leftf = centerf - TWOPI / ndata;
	rightf = centerf + TWOPI / ndata;
	
	f[1] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);
	
	amph(&A[1], &psi[1], f[1], x, y, ndata);
	
	for(j=1;j<=ndata;j++){
	  xdata[j] -= A[1]*cos( f[1]*(j-1) + psi[1] );
	  ydata[j] -= A[1]*sin( f[1]*(j-1) + psi[1] );
	}

	window(x, y, xdata, ydata, ndata);

	power(powsd, x, y, ndata); 
      }   

    else 
      centerf = freq[1][1];

    leftf = centerf - TWOPI / ndata;
    rightf = centerf + TWOPI / ndata;

    /* DETERMINE THE FIRST FREQUENCY */
    f[1] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);
    
    /* COMPUTE AMPLITUDE AND PHASE */
    amph(&A[1], &psi[1], f[1], x, y, ndata);
    
    /* SUBSTRACT THE FIRST HARMONIC FROM THE SIGNAL */
    for(j=1;j<=ndata;j++){
      xdata[j] -= A[1]*cos( f[1]*(j-1) + psi[1] );
      ydata[j] -= A[1]*sin( f[1]*(j-1) + psi[1] );
    }    
    
    /* HERE STARTS THE MAIN LOOP  *************************************/ 
    
    Q[1][1] = 1;
    alpha[1][1] = 1;
    
    for(m=2;m<=nfreq;m++){
      
      /* MULTIPLY SIGNAL BY WINDOW FUNCTION */
      window(x, y, xdata, ydata, ndata);
      
      /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
      power(powsd, x, y, ndata);
      
      if(l==1){
	
	centerf = bracket(powsd, ndata);

	leftf = centerf - TWOPI / ndata;
	rightf = centerf + TWOPI / ndata;

	f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);

	/* CHECK WHETHER THE NEW FREQUENCY IS NOT TOO CLOSE TO ANY PREVIOUSLY
	   DETERMINED ONE */
	nearfreqflag = 0.;
	for(k=1;k<=m-1;k++)
	  if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata )   nearfreqflag = 1; 
	    
	/* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
	while(f[m] < minfreq || f[m] > maxfreq || nearfreqflag == 1){
	  
	  /* IF NO, SUBSTRACT IT FROM THE SIGNAL */
	  leftf = centerf - TWOPI / ndata;
	  rightf = centerf + TWOPI / ndata;
	  
	  f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);
	  
	  amph(&A[m], &psi[m], f[m], x, y, ndata);
	  
	  for(j=1;j<=ndata;j++){
	    xdata[j] -= A[m]*cos( f[m]*(j-1) + psi[m] );
	    ydata[j] -= A[m]*sin( f[m]*(j-1) + psi[m] );
	  }
	  
	  /* AND RECOMPUTE THE NEW ONE */
	  window(x, y, xdata, ydata, ndata);
	  
	  power(powsd, x, y, ndata); 
	  
	  centerf = bracket(powsd, ndata); 

	  leftf = centerf - TWOPI / ndata;
	  rightf = centerf + TWOPI / ndata;
	  
	  f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);
	  
	  nearfreqflag = 0.;
	  for(k=1;k<=m-1;k++)
	    if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata )   nearfreqflag = 1; 

	}   

      } else {  
	
	centerf = freq[1][m];
	
	leftf = centerf - TWOPI / ndata;
	rightf = centerf + TWOPI / ndata;
	
	/* DETERMINE THE NEXT FREQUENCY */
	f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);
	
      }

      /* COMPUTE ITS AMPLITUDE AND PHASE */
      amph(&A[m], &psi[m], f[m], x, y, ndata);
      
      
      /* EQUATION (3) in Sidlichovsky and Nesvorny (1997) */
      Q[m][m] = 1;
      for(j=1;j<=m-1;j++){
	fac = (f[m] - f[j]) * (ndata - 1.) / 2.;
	Q[m][j] = sin(fac)/fac * PI*PI / (PI*PI - fac*fac);
	Q[j][m] = Q[m][j];
      }
      
      /* EQUATION (17) */
      for(k=1;k<=m-1;k++){
	B[k] = 0;
	for(j=1;j<=k;j++)
	  B[k] += -alpha[k][j]*Q[m][j];
      }
      
      /* EQUATION (18) */
      alpha[m][m] = 1;
      for(j=1;j<=m-1;j++)
	alpha[m][m] -= B[j]*B[j];
      alpha[m][m] = 1. / sqrt(alpha[m][m]);
      
      
      /* EQUATION (19) */
      for(k=1;k<=m-1;k++){
	alpha[m][k] = 0;
	for(j=k;j<=m-1;j++)
	  alpha[m][k] += B[j]*alpha[j][k];
	alpha[m][k] = alpha[m][m]*alpha[m][k];
      }
      
      /* EQUATION (22) */
      for(i=1;i<=ndata;i++){
	xsum=0; ysum=0;
	for(j=1;j<=m;j++){
	  fac = f[j]*(i-1) + (f[m]-f[j])*(ndata-1.)/2. + psi[m];
	  xsum += alpha[m][j]*cos(fac);
	  ysum += alpha[m][j]*sin(fac);
	}
	xdata[i] -= alpha[m][m]*A[m]*xsum;
	ydata[i] -= alpha[m][m]*A[m]*ysum;
      }
    }
    
    /* EQUATION (26) */
    for(k=1;k<=nfreq;k++){
      xsum=0; ysum=0;
      for(j=k;j<=nfreq;j++){
	fac = (f[j]-f[k])*(ndata-1.)/2. + psi[j];
	xsum += alpha[j][j]*alpha[j][k]*A[j]*cos(fac);
	ysum += alpha[j][j]*alpha[j][k]*A[j]*sin(fac);
      }
      A[k] = sqrt(xsum*xsum + ysum*ysum);
      psi[k] = atan2(ysum,xsum);
    }
    
    /* REMEMBER THE COMPUTED VALUES FOR THE FMFT */
    for(k=1;k<=nfreq;k++){
      freq[l][k] = f[k];
      amp[l][k] = A[k];
      phase[l][k] = psi[k];
    }
  }

 

  /* RETURN THE FINAL FREQUENCIES, AMPLITUDES AND PHASES */ 

  
   for(k=1;k<=nfreq;k++){
     signal1[k-1].freq = freq[1][k];            
     signal1[k-1].amp = amp[1][k];
     signal1[k-1].phase = phase[1][k];
 
     if(signal1[k-1].phase < -PI) signal1[k-1].phase += TWOPI;
     if(signal1[k-1].phase >= PI) signal1[k-1].phase -= TWOPI;
   }
   
   if(flag==2 || flag==3){
 
 
     for(k=1;k<=nfreq;k++){
       signal2[k-1].freq = freq[1][k] + (freq[1][k] - freq[2][k]);            
       signal2[k-1].amp = amp[1][k] + (amp[1][k] - amp[2][k]);
       signal2[k-1].phase = phase[1][k] + (phase[1][k] - phase[2][k]);
       
       if(signal2[k-1].phase < -PI) signal2[k-1].phase += TWOPI;
       if(signal2[k-1].phase >= PI) signal2[k-1].phase -= TWOPI;
     }
   
   if(flag==3){
     for(k=1;k<=nfreq;k++){
       
       signal3[k-1].amp = freq[1][k];
       if(fabs((fac = freq[2][k] - freq[3][k])/freq[2][k]) > FMFT_TOL)
 	signal3[k-1].freq += DSQR(freq[1][k] - freq[2][k]) / fac;
       else 
 	signal3[k-1].freq += freq[1][k] - freq[2][k]; 
 
       signal3[k-1].amp = amp[1][k];
       if(fabs((fac = amp[2][k] - amp[3][k])/amp[2][k]) > FMFT_TOL)
 	signal3[k-1].amp += DSQR(amp[1][k] - amp[2][k]) / fac;
       else
 	signal3[k-1].amp += amp[1][k] - amp[2][k]; 
 
       signal3[k].phase = phase[1][k];
       if(fabs((fac = phase[2][k] - phase[3][k])/phase[2][k]) > FMFT_TOL)
 	signal3[k-1].phase += DSQR(phase[1][k] - phase[2][k]) / fac;
       else
 	signal3[k-1].phase += phase[1][k] - phase[2][k]; 
 
       if(signal3[k-1].phase < -PI) signal3[k-1].phase += TWOPI;
       if(signal3[k-1].phase >= PI) signal3[k-1].phase -= TWOPI;
     }
   }
   }
 
   /* SORT THE FREQUENCIES IN DECREASING ORDER OF AMPLITUDE */
 
   int cmpfunc (const void * a, const void * b){
    return ( (*(struct component *)b).amp > (*(struct component *)a).amp);
   }
 
/*   if(flag==1) */
     qsort(signal1,nfreq,sizeof(struct component), cmpfunc);
   
/*   if(flag > 1){*/
     qsort(signal2,nfreq,sizeof(struct component), cmpfunc);
/*   }*/
 
/*   if(flag==3){*/
     qsort(signal3,nfreq,sizeof(struct component), cmpfunc);
/*   }*/
  /* FREE THE ALLOCATED VARIABLES */
/*  free_dvector(xdata, 1, ndata);*/
/*  free_dvector(ydata, 1, ndata);*/

  free_dvector(x, 1, ndata);
  free_dvector(y, 1, ndata);
  free_vector(powsd, 1, ndata);
  
  free_dmatrix(freq, 1, 3*flag, 1, nfreq); 
  free_dmatrix(amp, 1, 3*flag, 1, nfreq);
  free_dmatrix(phase, 1, 3*flag, 1, nfreq);

  free_dvector(f, 1, nfreq);
  free_dvector(A, 1, nfreq);
  free_dvector(psi, 1, nfreq);
 
  free_dmatrix(Q, 1, nfreq, 1, nfreq); 
  free_dmatrix(alpha, 1, nfreq, 1, nfreq);
  free_dvector(B, 1, nfreq);

  return 1;
}


void window(double *x, double *y, double *xdata, double *ydata, size_t ndata)

/* MULTIPLIES DATA BY A WINDOW FUNCTION */      
{  
  long j;
  double window;
  
  for(j=1;j<=ndata;j++) {

    window = TWOPI*(j-1) / (ndata-1);
    window = (1. - cos(window)) / 2.;

    x[j] = xdata[j]*window;
    y[j] = ydata[j]*window;

  }
}


void power(float *powsd, double *x, double *y, size_t ndata)

/* REARRANGES DATA FOR THE FAST FOURIER TRANSFORM, 
CALLS FFT AND RETURNS POWER SPECTRAL DENSITY */

/* todo: use the gsl_fft_complex_backward */
  /* require gsl_ff_complex_watable and complex_workspac
     but trival to allocate and deallocate */

{
  long j;
  double z[2*ndata];


   /* I use a simple "double" and therefore does not need to be allocated
    * or de-allocated, but beware that indexing begins at zero */

  for(j=1;j<=ndata;j++){
    z[2*j-2] = x[j];
    z[2*j-1] = y[j];
  }
 
  gsl_fft_complex_radix2_backward (z,1,ndata);
 
  for(j=1;j<=ndata;j++)
    powsd[j] = FSQR(z[2*j-2]) + FSQR(z[2*j-1]);

}

double bracket(float *powsd, size_t ndata)

/* mc: as far as I can tell this does not come from NR */
/* FINDS THE MAXIMUM OF THE POWER SPECTRAL DENSITY  */ 

{
  long j, maxj;
  double freq, maxpow;

  maxj = 0;
  maxpow = 0;
 
  for(j=2;j<=ndata/2-2;j++)
    if(powsd[j] > powsd[j-1] && powsd[j] > powsd[j+1])
      if(powsd[j] > maxpow){ 
	maxj = j;
        maxpow = powsd[j];
      }  

  for(j=ndata/2+2;j<=ndata-1;j++)
    if(powsd[j] > powsd[j-1] && powsd[j] > powsd[j+1])
      if(powsd[j] > maxpow){ 
	maxj = j;
        maxpow = powsd[j];
      }  

  if(powsd[1] > powsd[2] && powsd[1] > powsd[ndata])
    if(powsd[1] > maxpow){ 
      maxj = 1;
      maxpow = powsd[1];
    }

  if(maxpow == 0) nrerror("DFT has no maximum ...");

  if(maxj < ndata/2) freq = -(maxj-1);  
  if(maxj > ndata/2) freq = -(maxj-ndata-1);

  return (TWOPI*freq / ndata);

  /* negative signs and TWOPI compensate for the Numerical Recipes 
     definition of the DFT */
}
    
#define GOLD_R 0.61803399
#define GOLD_C (1.0 - GOLD_R)

double golden(double (*f)(double, double *, double *, size_t), 
	      double ax, double bx, double cx,
	      double xdata[], double ydata[], size_t n)

     /* calculates the maximum of a function bracketed by ax, bx and cx */

{
  double f1,f2,x0,x1,x2,x3;

  x0=ax;
  x3=cx;

  if(fabs(cx-bx) > fabs(bx-ax)){
    x1 = bx;
    x2 = bx + GOLD_C*(cx-bx);
  } else {
    x2 = bx;
    x1 = bx - GOLD_C*(bx-ax);
  }

  f1 = (*f)(x1, xdata, ydata, n);
  f2 = (*f)(x2, xdata, ydata, n);

  while(fabs(x3-x0) > FMFT_TOL*(fabs(x1)+fabs(x2))){
    if(f2 > f1){
      SHFT4(x0,x1,x2,GOLD_R*x1+GOLD_C*x3);
      SHFT3(f1,f2,(*f)(x2, xdata, ydata, n));
    } else {
      SHFT4(x3,x2,x1,GOLD_R*x2+GOLD_C*x0);
      SHFT3(f2,f1,(*f)(x1, xdata, ydata, n));
    }
  }

  if(f1>f2) return x1;
  else return x2;
}

void amph(double *amp, double *phase, double freq, 
	  double xdata[], double ydata[], size_t ndata){

  /* CALCULATES THE AMPLITUDE AND PHASE */

  double xphi,yphi;
  
  xphi = yphi = 0;
  
  phifun(&xphi, &yphi, freq, xdata, ydata, ndata);
  
  *amp = sqrt(xphi*xphi + yphi*yphi);
  *phase = atan2(yphi, xphi);
}

double phisqr(double freq, double xdata[], double ydata[], size_t ndata)

/* COMPUTES A SQUARE POWER OF THE FUNCTION PHI */

{	
  double xphi,yphi;
  
  xphi = yphi = 0;
  
  phifun(&xphi, &yphi, freq, xdata, ydata, ndata);
  
  return xphi*xphi + yphi*yphi;
}

void phifun(double *xphi, double *yphi, double freq,  
	      double xdata[], double ydata[], long n)

     /* COMPUTES THE FUNCTION PHI */   

{
  long i, j, nn;
  double c, s, *xdata2, *ydata2;
  
  xdata2 = dvector(1, n);
  ydata2 = dvector(1, n);
  
  xdata2[1] = xdata[1] / 2; ydata2[1] = ydata[1] / 2;
  xdata2[n] = xdata[n] / 2; ydata2[n] = ydata[n] / 2;

  for(i=2;i<=n-1;i++){
    xdata2[i] = xdata[i];
    ydata2[i] = ydata[i];
  }

  nn = n;

  while(nn != 1){
    
    nn = nn / 2;
    
    c = cos(-nn*freq);
    s = sin(-nn*freq);
   
    for(i=1;i<=nn;i++){
      j=i+nn;
      xdata2[i] += c*xdata2[j] - s*ydata2[j];
      ydata2[i] += c*ydata2[j] + s*xdata2[j];
    }

  }
  
  *xphi = 2*xdata2[1] / (n-1);
  *yphi = 2*ydata2[1] / (n-1);

  free_dvector(xdata2,1,n);
  free_dvector(ydata2,1,n);

}

#define SORT_M 7 
#define SORT_NSTACK 50

