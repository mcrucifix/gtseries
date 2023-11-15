# I need to add lots of comments about licencing
# the fact that the hanning window is now made here
# check that it provides the same output as the original code
# espacially with the new window

#' Modified Fourier transform (Real C)
#'
#' Implementation of the the Frequency Modified Fourier Transform 
#' (Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137). 
#' Given a quasi--periodic complex signal X + iY, the algorithm 
#' estimates the frequencies (f_j), amplitudes (A_j) and phases 
#' (psi_j) in its decomposition. 
#' @useDynLib gtseries
#' @importFrom stats start
#' @importFrom stats as.ts
#' @param xdata The data provided either as a time series (advised), or as a vector. 
#' may be complex
#' @param minfreq,maxfreq If provided, bracket the frequencies to be probed. Note this are
#'        more precisely angular velocities (2\pi / period), expressed in time-inverse units
#'        with the time resolution encoded in `xdata` if the latter is a time series. 
#' @param flag :
#'      Modified Fourier Transform                  if   flag = 1;
#'      Frequency Modified Fourier Transform        if   flag = 2;
#'      FMFT with additional non-linear correction  if   flag = 3
#' (while the first algorithm is app. 3 times faster than the third one, 
#' the third algorithm should be in general much more precise).  
#' The computed frequencies are in the range given by minfreq and maxfreq.
#' @param nfreq is the number of frequencies returned, must be smaller that the length of  xdata.
#' @return STILL NEED TO BE DESCRIBED
#' @author Michel Crucifix for the R code, and David Nesvorny for most of the supporting C code doing the
#' actual computations
#' \insertRef{sidlichovsky97aa}{gtseries}
#
#' @export mfft_real_C
mfft_real_C <- function(xdata, minfreq=NULL, maxfreq=NULL, flag=1, nfreq=30)
{

  xdata = stats::as.ts(xdata)
  dt = deltat(xdata)
  startx = stats::start(xdata)[1]
  ndata = length(xdata);


  if (!(flag == 1 || flag == 2 || flag == 3)) stop("flag must be either 1, 2, or 3")

#   if (window=="hanning")
#   {
#     hanning = (1-cos(2*pi*(seq(ndata)-1)/(ndata-1)))/2.
#     windowing
#     xdata = xdata*hanning
#     ydata = ydata*hanning
#   } else if (window=="none") {} else stop ("this window is not supported, sorry")
#     padding
#     N = 2^ceiling(log(ndata)/log(2))
#     if (N > ndata){
#       xdata[(ndata+1):N] = 0
#       ydata[(ndata+1):N] = 0
#     }
#     

    if (is.null(minfreq)){
      my_minfreq = -0.001
      my_maxfreq = pi} else {
      my_minfreq = minfreq * dt
      my_maxfreq = maxfreq * dt}
 
     signal1 = as.double(rep(0,3*nfreq))
     signal2 = as.double(rep(0,3*nfreq))
     signal3 = as.double(rep(0,3*nfreq))

     Infreq=as.integer(nfreq)
     Iminfreq = as.double(my_minfreq)
     Imaxfreq = as.double(my_maxfreq)
     Ixdata = as.double(xdata)

     OUT = .C("fmft_real", Infreq, Iminfreq, Imaxfreq, as.integer(flag), 
              as.integer(ndata), Ixdata, signal1,signal2,signal3, DUP=TRUE)

     OUTVEC <- t(matrix(OUT[[7+flag]], 3, nfreq))

     Freq <- OUTVEC[seq(nfreq) ]/dt
     Ampl <- OUTVEC[nfreq+seq(nfreq) ] 
     Phase <- OUTVEC[2*nfreq+seq(nfreq) ] - startx*Freq

     # if this is a real vector, there will be positive and negative frequencies corresponding to the same amplitude
     # take care and verify that this actually works this way

     return(data.frame(Freq=Freq, Ampl=Ampl, Phase=Phase))

  }
#

