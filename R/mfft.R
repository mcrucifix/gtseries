# I need to add lots of comments about licencing
# the fact that the hanning window is now made here
# check that it provides the same output as the original code
# espacially with the new window

#' Modified Fourier transform 
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
#' @export mfft
mfft <- function(xdata, minfreq=NULL, maxfreq=NULL, flag=1, nfreq=30)
{

  xdata = stats::as.ts(xdata)
  dt = deltat(xdata)
     print('deltat')
     print(dt)
  startx = stats::start(xdata) 
  ydata = Im(xdata)
  xdata = Re(xdata)
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
      my_minfreq = -pi
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
     Iydata = as.double(ydata)

     OUT = .C("fmft", Infreq, Iminfreq, Imaxfreq, as.integer(flag), 
              as.integer(ndata), Ixdata, Iydata, signal1,signal2,signal3, DUP=TRUE)

     OUTVEC <- t(matrix(OUT[[7+flag]], 3, nfreq))
  

     Freq <- OUTVEC[seq(nfreq) ]/dt
     Ampl <- OUTVEC[nfreq+seq(nfreq) ] 
     Phase <- OUTVEC[2*nfreq+seq(nfreq) ] - startx*Freq

     # if this is a real vector, there will be positive and negative frequencies corresponding to the same amplitude
     # take care and verify that this actually works this way

     return(data.frame(Freq=Freq, Ampl=Ampl, Phase=Phase))

  }
#

#
#                         return(OUTVEC)
# dyn.load('fmft.so')
# ndata=8192

#for (exp in c('d')){
#  namein <- file.path('La2010',sprintf("La2010%s_alkhqp3L.dat",exp));
#  nameout <- sprintf("La2010%s_out.RData",exp);
#  A <- read.table(namein)
#  deltat=1000
#  colnames(A) <- c('times','a','l','k','h','q','p')
#  # ndata = as.integer64(ndata)
#
#  times <- seq(ndata)
#
#  nfreq <- as.integer(30)
#  flag =  2 
#
#  minfreq <- as.double ( -1e6 / 180. / 3600. * pi /5. ) ;
#  maxfreq <- as.double ( 1.e6 / 180. / 3600. * pi /5. ) ;
#
#  signal1 = as.double (rep(0, nfreq*3))
#  signal2 = as.double (rep(0, nfreq*3))
#  signal3 = as.double (rep(0, nfreq*3))
#
#  ntrunks <- floor(length(A$times)/ndata)
#
#  # EMFFT <- sapply (seq(0,floor(length(A$times)/ndata)-1), function(i)
#
#
#  timemin <- sapply (seq(0,ntrunks-1), function(i) { A$times [ndata* i + 1] });
#  timemax <- sapply (seq(0,ntrunks-1), function(i) { A$times [ndata* i + ndata ] });
#
#  EMFFT <- sapply (seq(0,ntrunks-1), function(i)
#                       {
#                         xdata <- A$h[ndata*i + (1:ndata)]
#                         ydata <- A$k[ndata*i + (1:ndata)]
#
#                         OUT = .C("fmft", nfreq, minfreq, maxfreq, as.integer(flag), as.integer(ndata), xdata, ydata, signal1,signal2,signal3, DUP=FALSE)
#
#                         # OUTVEC <- matrix(OUT[[8]], 10,3)
#                         OUTVEC <- t(matrix(OUT[[8]], 3, nfreq))
#
#                         return(OUTVEC)
#
#                       })
#
#
#  # repack
#
#  # repack
#
#
#  Freq <- EMFFT[seq(nfreq), ]*360*60*60/2/pi/deltat
#  Ampl <- EMFFT[nfreq+seq(nfreq), ]
#  Phas <- EMFFT[2*nfreq+seq(nfreq), ]
#
#
#  snake <- function(x1, x2, tol=2.e-1, mynfreq=nfreq){
#    xnew <- c();
#    xamp <- c();
#    xpha <- c();
#    xf1  <- x1[seq(mynfreq)]
#    xf2  <- x2[seq(mynfreq)]
#    xa1  <- x1[mynfreq+seq(mynfreq)]
#    xa2  <- x2[mynfreq+seq(mynfreq)]
#    xp1  <- x1[2*mynfreq+seq(mynfreq)]
#    xp2  <- x2[2*mynfreq+seq(mynfreq)]
#
#    while(length(xf1) > 0) {
#      i <- which.min(abs(xf2 - xf1[1])+0*abs(xa2-xa1[1]));
#      xnew <- c(xnew, xf2[i])
#      xamp <- c(xamp, xa2[i])
#      xpha <- c(xpha, xp2[i])
#      xf2 <- xf2[-i]
#      xf1 <- xf1[-1]
#      xa2 <- xa2[-i]
#      xa1 <- xa1[-1]
#      xp2 <- xp2[-i]
#      xp1 <- xp1[-1]
#    }
#    return(c(xnew, xamp, xpha))
#  }
#
#  mylist <- apply(EMFFT,2,as.numeric, simplify=FALSE)
#
#  terms <- simplify2array(Reduce(snake, mylist, mylist[[1]], accumulate=TRUE))
#
#  myfreq <- terms[seq(nfreq),-1] * 360*60*60/2/pi/deltat
#  myamp <- terms[nfreq+seq(nfreq),-1]
#  myphase <- terms[2*nfreq+seq(nfreq),-1]
#
#
#  rownames (myfreq) <- timemin
#  rownames (myamp)  <- timemin
#  rownames (myphase)<- timemin
#
#
#  OUT = list(Freq=myfreq, Amp=myamp, Phase=myphase)
#  save(OUT, file=nameout)
#}:
