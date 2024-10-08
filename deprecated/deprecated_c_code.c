int fmft_real(int *localnfreq, double *localminfreq, double *localmaxfreq, int *localflag, 
	 int *localndata, double *localxdata, 
   struct component *signal1, struct component *signal2, struct component *signal3);


int fmft_real(int *localnfreq, double *localminfreq, double *localmaxfreq, int *localflag, 
	 int *localndata, double *localxdata,
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
  double *xdata,  *x, *y;
  double centerf, leftf, rightf, facplus, facminus, fac,
         sinplus, cosplus, sinminus, cosminus, factemp,
         xsum, ysum;
  double **freq, **amp, **phase, *f, *A, *Ac, *As, *psi;
  double **Q, **alpha, *B;

  FILE *fp;

  int nfreq = *localnfreq;
  double minfreq = *localminfreq;
  double maxfreq = *localmaxfreq;
  int flag = *localflag;
  size_t ndata = *localndata;
  size_t ndata_real = (ndata+1)/2 ; 
  fastflag = isPowerofTwo(ndata) ;

  if (fastflag)
    (printf("ndata is power of two: we will be faster ! \n"));
  if (ndata <= 2){
    printf("at least 2 data needed - output non-reliable"); return(0);
  }
  if (ndata <= nfreq){
    printf("nfreq must be smaller than nata"); return(0);
  }

/*  printf("prelimarg %d, %zu %d", nfreq, ndata, flag);*/

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
  As = dvector(1, nfreq);
  Ac = dvector(1, nfreq);
  psi = dvector(1, nfreq);

  
  Q = dmatrix(1, 2*nfreq, 1, 2*nfreq); 
  alpha = dmatrix(1, 2*nfreq, 1, 2*nfreq);
  B = dvector(1, 2*nfreq);
 

  /* 1 LOOP FOR MFT, 2 LOOPS FOR FMFT, 3 LOOPS FOR NON-LINEAR FMFT */

  for(l=1; l<=flag; l++){
 
    if(l==1){
      xdata = localxdata -1;  // -1 because dvector vs *double
/*      ydata = localydata -1;*/
      /* SEPARATE REAL AND IMAGINERY PARTS */ 
/*      for(j=1;j<=ndata;j++){*/
/*  xdata[j] = localxdata[j-1];*/
/*  ydata[j] = localydata[j-1];*/

    } else {

       /* GENERATE THE QUASIPERIODIC FUNCTION COMPUTED BY MFT */
      for(i=1;i<=ndata;i++){
	xdata[i] = 0; 
	for(k=1;k<=nfreq;k++){
	  xdata[i] += amp[l-1][k]*cos(freq[l-1][k]*(i-1) + phase[l-1][k]);
/*    ydata[i] += amp[l-1][k]*sin(freq[l-1][k]*(i-1) + phase[l-1][k]);*/
	}
      }

    }
  

    /* MULTIPLY THE SIGNAL BY A WINDOW FUNCTION, STORE RESULT IN x AND y */
    window_real(x, xdata, ndata);
    
    /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
    power_real(powsd, x, ndata);


    if(l==1)  {

	printf("l=1 ; start the while loop \n");
      /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
      while((centerf = bracket_real(powsd, ndata)) < minfreq || centerf > maxfreq) {

	printf("centerf = %.2f \n",centerf);
	/* IF NO, SUBSTRACT IT FROM THE SIGNAL */
	leftf = centerf - TWOPI / ndata;
	rightf = centerf + TWOPI / ndata;
	
	f[1] = golden_real(phisqr_real, leftf, centerf, rightf, x, ndata);
	
	printf("f[1] = %.2f \n",centerf);
	amph_real(&A[1], &Ac[1], &As[1], &psi[1], f[1], x, ndata);
	printf("&A[1] = %.2f ; Ac = %.2f ; As = %.2f\n",&A[1], &Ac[1], &As[1]);
	printf("------\n",centerf);
	
	for(j=1;j<=ndata;j++){
/*    xdata[j] -= A[1]*cos( f[1]*(j-1) + psi[1] );*/
	  xdata[j] -= Ac[1]*cos( f[1]*(j-1)  );
	  xdata[j] -= As[1]*sin( f[1]*(j-1)  );
	}


	window_real(x, xdata, ndata);
	power_real(powsd, x, ndata); 
      }   }

    else 
      centerf = freq[1][1];

    leftf = centerf - TWOPI / ndata;
    rightf = centerf + TWOPI / ndata;

    /* DETERMINE THE FIRST FREQUENCY */
    f[1] = golden_real(phisqr_real, leftf, centerf, rightf, x, ndata);
    
    /* COMPUTE AMPLITUDE AND PHASE */
    amph_real(&A[1], &Ac[1], &As[1], &psi[1], f[1], x, ndata);
    
    /* SUBSTRACT THE FIRST HARMONIC FROM THE SIGNAL */
    for(j=1;j<=ndata;j++){
      xdata[j] -= Ac[1]*cos( f[1]*(j-1)  );
      xdata[j] -= As[1]*sin( f[1]*(j-1)  );
    }    
    /* HERE STARTS THE MAIN LOOP  *************************************/ 
    
    printf("start the main loop \n");
    Q[1][1] = 1;
    alpha[1][1] = 1;
    
    for(m=2;m<=nfreq;m++){
    printf("m = %i \n", m );
      /* MULTIPLY SIGNAL BY WINDOW FUNCTION */
      window_real(x, xdata, ndata);
      
      /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
      power_real(powsd, x, ndata);
      
      if(l==1){
	
	centerf = bracket_real(powsd, ndata);

	leftf = centerf - TWOPI / ndata;
	rightf = centerf + TWOPI / ndata;

	f[m] = golden_real(phisqr_real, leftf, centerf, rightf, x, ndata);

	/* CHECK WHETHER THE NEW FREQUENCY IS NOT TOO CLOSE TO ANY PREVIOUSLY
	   DETERMINED ONE */
	nearfreqflag = 0;
	for(k=1;k<=m-1;k++) if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata )   nearfreqflag = k; 
	    
	/* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
	while(f[m] < minfreq || f[m] > maxfreq || nearfreqflag > 0){
	  
	  printf("centerf = %.2f \n",centerf);
	  /* IF NO, SUBSTRACT IT FROM THE SIGNAL */
	  leftf = centerf - TWOPI / ndata;
	  rightf = centerf + TWOPI / ndata;
	  
	  f[m] = golden_real(phisqr_real, leftf, centerf, rightf, x, ndata);
	  
	printf("f[%i] = %.2f (minfreq= %.2f, maxfreq=%.2f) \n",m, f[m], minfreq, maxfreq);
	  amph_real(&A[m], &Ac[m], &As[m], &psi[m], f[m], x,  ndata);
	printf("&A[1] = %.2f ; Ac = %.2f ; As = %.2f\n",&A[1], &Ac[1], &As[1]);
	printf("------\n",centerf);
	  
	  for(j=1;j<=ndata;j++){
	    xdata[j] -= Ac[m]*cos( f[m]*(j-1)  );
	    xdata[j] -= As[m]*sin( f[m]*(j-1)  );
	  }
	  
	  /* AND RECOMPUTE THE NEW ONE */
	  window_real(x, xdata, ndata);
	  
	  power_real(powsd, x, ndata); 
	  
	  centerf = bracket_real(powsd, ndata); 

	  leftf = centerf - TWOPI / ndata;
	  rightf = centerf + TWOPI / ndata;
	  
	  f[m] = golden_real(phisqr_real, leftf, centerf, rightf, x, ndata);
	  
	  nearfreqflag = 0.;
	  for(k=1;k<=m-1;k++)
	    if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata )   nearfreqflag = 1; 

	}   

      } else {  
	
	centerf = freq[1][m];
	
	leftf = centerf - TWOPI / ndata;
	rightf = centerf + TWOPI / ndata;
	
	/* DETERMINE THE NEXT FREQUENCY */
	f[m] = golden_real(phisqr_real, leftf, centerf, rightf, x,  ndata);
	
      }

      /* COMPUTE ITS AMPLITUDE AND PHASE */
      amph_real(&A[m], &Ac[m], &As[m],  &psi[m], f[m], x, ndata);
      
      /* EQUATION (3) in Sidlichovsky and Nesvorny (1997) */
	    facplus = (f[m] + f[m]) * (ndata - 1.) / 2.;
	    facminus = 0.;

      facplus= sin(facplus)/facplus * PI*PI / (PI*PI - facplus*facplus);
      facminus= 1.;

      sinplus = 0.5 * sin(2*f[m]) * facplus;
      sinminus = 0.;
      cosplus = 0.5 * cos(2*f[m]) * facplus;
      cosminus =  facplus;
    	Q[2*m-1][2*m-1] = (cosminus + cosplus);
    	Q[2*m-1][2*m]   = sinplus;
    	Q[2*m][2*m-1]   = sinplus;
    	Q[2*m][2*m]     = (cosminus - cosplus);


      for(j=1;j<=m-1;j++){

	facplus = (f[m] + f[j]) * (ndata - 1.) / 2.;
	facminus = (f[m] - f[j]) * (ndata - 1.) / 2.;
  facminus= sin(facminus)/facminus * PI*PI / (PI*PI - facminus*facminus);
  facplus= sin(facplus)/facplus * PI*PI / (PI*PI - facplus*facplus);


  sinplus  = 0.5 * sin(f[m] + f[j]) * facplus ; 
  sinminus = 0.5 * sin(f[m] - f[j]) * facminus ; 
  cosplus  = 0.5 * sin(f[m] + f[j]) * facplus ; 
  cosminus = 0.5 * sin(f[m] + f[j]) * facminus ; 

  // coscos, cossin, sincos, sinsin
	Q[2*m-1][2*j-1] = (cosminus + cosplus);
	Q[2*m-1][2*j] =   (sinplus - sinminus);
	Q[2*m][2*j-1] =   (sinplus + sinminus);
	Q[2*m][2*j] =     (cosminus - cosplus);

  // the symmetrics are, indeed, symmetric
  //
	Q[2*m-1][2*j-1] = Q[2*j-1][2*m-1];
	Q[2*m-1][2*j  ] = Q[2*j-1][2*m  ];
	Q[2*m  ][2*j-1] = Q[2*j  ][2*m-1];
	Q[2*m  ][2*j  ] = Q[2*j  ][2*m  ];

      }
      
      /* EQUATION (17) */
      for(k=1;k<=2*m-1;k++){
	B[k] = 0;
	for(j=1;j<=k;j++)
	  B[k] += -alpha[k][j]*Q[m][j];
      }

      /* EQUATION (18) */
      alpha[m][m] = 1;
      for(j=1;j<=2*m-1;j++)
	alpha[m][m] -= B[j]*B[j];
      alpha[m][m] = 1. / sqrt(alpha[m][m]);
      
      
      /* EQUATION (19) */
      for(k=1;k<=2*m-1;k++){
	alpha[m][k] = 0;
	for(j=k;j<=2*m-1;j++)
	  alpha[m][k] += B[j]*alpha[j][k];
	alpha[m][k] = alpha[m][m]*alpha[m][k];
      }
 
/* ICICICICICI */ 

      /* EQUATION (22) */
      for(i=1;i<=ndata;i++){
	xsum=0; ysum=0;
  /* on est a l'equation 21
   * la difficulte est que f_m est en fait f_(2m)
   * il faut retnancher deux composantet
   * f[2m] = f_[2m-1] - a(2m)(2m) f_[2m-1]*sum1
   * f_[2m-1] = f_[2m-2] - a(2m-1)(2m-1) f_[2m-2]*sum2 */
	for(j=1;j<=2*m;j=j+2){
    fac = f[j]*(i-1) ;
	  xsum += alpha[2*m-1][2*j-1]*cos(fac);
	  xsum += alpha[2*m-1][2*j]*sin(fac);
	  ysum += alpha[2*m][2*j-1]*cos(fac);
	  ysum += alpha[2*m][2*j]*sin(fac);
	}
	xdata[i] -= alpha[2*m-1][2*m-1]*Ac[m]*xsum;
	xdata[i] -= alpha[2*m][2*m]*As[m]*ysum;
	xdata[i] -= alpha[2*m][2*m-1]*As[m]*sin(fac);
      }
    }
    
    /* EQUATION (26) */
    for(k=1;k<=nfreq;k++){
      xsum=0; ysum=0;
      for(j=k;j<=nfreq;j++){
      	xsum += alpha[2*j-1][2*j-1]*alpha[2*j-1][2*k-1]*Ac[j];
      	xsum += alpha[2*j-1][2*j-1]*alpha[2*j-1][2*k]*Ac[j];
      	ysum += alpha[2*j][2*j]*alpha[2*j][2*k-1]*As[j];
      	ysum += alpha[2*j][2*j]*alpha[2*j][2*k]*As[j];
      }
        xsum += alpha[2*k-1][2*k-1]*alpha[2*k-1][2*k]*Ac[k];
       A[k] = sqrt(xsum*xsum + ysum*ysum);
       Ac[k] = xsum;
       As[k] = ysum;
       psi[k] = -atan2(ysum,xsum);
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
  free_dvector(Ac, 1, nfreq);
  free_dvector(As, 1, nfreq);
  free_dvector(psi, 1, nfreq);
 
  free_dmatrix(Q, 1, 2*nfreq, 1, 2*nfreq); 
  free_dmatrix(alpha, 1, 2*nfreq, 1, 2*nfreq);
  free_dvector(B, 1, 2*nfreq);

  return 1;
}
