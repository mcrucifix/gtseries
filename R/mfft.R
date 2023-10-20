# I need to add lots of comments about licencing
# the fact that the hanning window is now made here
# check that it provides the same output as the original code
# espacially with the new window



mfft <- function(xdata, ydata=NULL, minfreq=NULL, maxfreq=NULL, flag=2, window="hanning", nfreq=30)
{

  # please at the moent do not play with 'flag'
  ndata = length(xdata);
  if (is.null(ydata)) {
    ydata = rep(0, ndata)} else {
      if (length(ydata) != ndata) stop("I am confused: you provide a ydata that does not match the length of xdata")
  }

  if (window=="hanning")
  {
    hanning = (1-cos(2*pi*(seq(ndata)-1)/(ndata-1)))/2.
    # windowing
    xdata = xdata*hanning
    ydata = ydata*hanning
  } else if (window=="none") {} else stop ("this window is not supported, sorry")
    # padding
    N = 2^ceiling(log(ndata)/log(2))
    if (N > ndata){
      xdata[(ndata+1):N] = 0
      ydata[(ndata+1):N] = 0
    }
    
    xdata = ts(xdata)
    ydata = ts(ydata, start=start(xdata), deltat=deltat(xdata))
    dt = deltat(xdata)

    if (is.null(minfreq)){
      my_minfreq = pi/ndata
      my_maxfreq = pi} else {
      my_minfreq = minfreq * dt
      my_maxfreq = maxfreq * dt}
 
    # to do: separate minfreq from maxfreq

     signal1 = as.double(rep(0,3*nfreq))
     signal2 = as.double(rep(0,3*nfreq))
     signal3 = as.double(rep(0,3*nfreq))
 

     print('length signal1')
     print(length(as.double(signal1)))
     print(length(as.double(signal2)))
     print(length(as.double(signal3)))
     print(nfreq)
  
     print(c("xdata:",xdata[1],xdata[2], xdata[3]))
     print(c("ydata:",ydata[1],ydata[2], ydata[3]))
   
     Infreq=as.integer(nfreq)
     Iminfreq = as.double(my_minfreq)
     Imaxfreq = as.double(my_maxfreq)
     Ixdata = as.double(xdata)
     Iydata = as.double(ydata)


     print('length xdata')
     print(length(Ixdata))
     print(length(Iydata))
     print(length(signal1))
     print(N)

     print("N")
     print(N)
     print(Ixdata)
     print("Iydata")
     print(Iydata)
#      print (c("fmft", Infreq, Iminfreq, Imaxfreq, as.integer(flag), as.integer(N), Ixdata, Iydata))
     OUT = .C("fmft", Infreq, Iminfreq, Imaxfreq, as.integer(flag), as.integer(N), Ixdata, Iydata, signal1,signal2,signal3, DUP=TRUE)
#      OUT = .C("fmft", nfreq, minfreq, maxfreq, as.integer(flag), as.integer(ndata), xdata, ydata, signal1,signal2,signal3, DUP=FALSE)

     OUTVEC <- t(matrix(OUT[[8]], 3, nfreq))
     print("OUTVEC")
     print(OUTVEC)

     Freq <- OUTVEC[seq(nfreq) ]/dt
     Ampl <- OUTVEC[nfreq+seq(nfreq) ] * N/ndata
     Phase <- OUTVEC[2*nfreq+seq(nfreq) ]

     # if this is a real vector, there will be positive and negative frequencies corresponding to the same amplitude
     # take care and verify that this actually works this way

     return(list(Freq=Freq, Ampl=Ampl, Phase=Phase))

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
#}
