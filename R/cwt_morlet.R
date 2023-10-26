#This program performs the Continuous Wavelet Transform (CWT) 
#of the input time series y. 
#It plots the series in normalized form 
#and displays the modulus (amplitude) of the CWT 
#in the time-period space. 
#The period is expressed in unit of time.
#
#
#   are in the vector named "period" 

#reference:
#
#Mallat, S. 1998: A wavelet Tour of Signal Processing. 
#Academic Press, 577 pp.
#
## for ridge extraction : 
## Tchawitchia P., Wavelets, functions and operators, ch. 3 in 
## Wavelets : theory and applications, Erlenbacher et al. eds, Oxford University Press 1996




col_wavelet <- colorRampPalette(c('darkblue','lightblue','grey','white','red'))(100)

search_ridge <- function (A,first_guess,k0=5.6,tol=1e-6, itmax=20,B=seq(along=A),window=c(0,Inf))
{
  n <- length(A)
  Amp <- A
  Amp[] <- NA
  Ome <- Amp
  Scale <- Amp
  deltat <- deltat(A)
   scale <- first_guess
   for (b in B)
   {
     delta = Inf
     it <- 0
     while ((delta > tol) & (it <= itmax)) {
         w <- cwt_morlet(A,k0=k0,calcmask=FALSE,scale=scale,deriv=TRUE)
         p <- attr(w,"deriv")[b]
         expected_scale = k0/Re(p)
         delta <- scale - expected_scale 
         scale <- expected_scale 
         it <- it+1
         if (is.na(scale)) {it = itmax}
       }
       if (it < itmax &  scale > window[1] & scale < window[2]  ) { 
           Amp[b] = w[b]
           Ome[b] = p 
           Scale[b] = scale
           }
       else {     Amp[b] = NA
                  Ome[b] = NA 
                  Scale[b] = NA 
                  if (scale > window[2]) scale = window[2] 
                  if (scale < window[1]) scale = window[1] 
            }      
   } 

    
   attr(Amp,"period") = 2*pi*deltat/Ome
   attr(Amp,"scale") = Scale
   attr(Amp,"k0") = k0
   attr(Amp,"class") = "ridge"
   Amp
}

enhance_ridge <- function(A,R,confidence=0.95,inter=0.08)
{
 Scale = attr(R,"scale")
 Amp = R
 period=matrix(NA,nrow=length(A),ncol=3)
 k0 = attr(R,"k0")
 deltat = deltat(A)


 width=log(qnorm(0.5+confidence/2)*sqrt(2)/k0+1)

 for (b in seq(along=R))
  {
    print(b)
    if (!is.na(Scale[b]))
    {
      scale=exp(log(Scale[b])+seq(-width*1.03214,+width,inter))
      ## 1.03214 is an empirical factor to take into account the ridge asymmetry 
      w <- cwt_morlet(A,k0=k0,calcmask=FALSE,scale=scale,deriv=TRUE)
      Amp[b]=sum(w[b,])*k0/sqrt(2*pi)*inter
      period[b,] = 2*pi*deltat/k0*c(Scale[b],min(scale),max(scale))
    }
   }
    attr(Amp,"period") <- period
    Amp
}


reconstruct_morlet  <- function(W,scales=c(-Inf,Inf), periods=NULL)
{
 if (!(attr(W,"class")=="wavelet")) stop ("object is not a wavelet transform")
 if (!(attr(W,"wavelet")=="morlet")) stop ("object is not a MORLET  wavelet transform")
 a <- attr(W,"scale")
 j <- which(a > scales[1] & a < scales[2])
 # if period is set, then overrides scales
 if (!is.null(periods))
 {
   a <- attr(W,"period")
   j <- which(a > periods[1] & a < periods[2])
 }

 inter <- attr(W,"parameters")$inter
 k0 <- attr(W,"parameters")$k0

 T <- rowSums(W[,j])*log(2)/inter/(sqrt(2*pi))*k0
 T <- ts(T, start=attr(W,'time')[1], deltat = diff(attr(W,'time'))[1])
}

cross_morlet <- function(A, B, ...)
{
  CA <- cwt_morlet(A, ...)
  CB <- cwt_morlet(B, ...)
  (CA * Conj(CB)) / (Mod(CA) * Mod(CB) )
}

#' Continous Morlet Wavelet Transform 

#' @importFrom stats fft deltat lm qnorm toeplitz ts time
#' @importFrom graphics axis image lines par text
#' @export cwt_morlet
cwt_morlet <- function (A,inter=20,k0=5.6,amin=1,amax=Inf,calcmask=TRUE,scale=NA,deriv=FALSE)
{
   y <- A
   xx <- stats::time(A)
   deltat<-deltat(A)
   ny<-length(y);

  if (is.na(scale)) {
    local({
      ny2<-round(ny/2);
      exp1<-log2(amin)+2;
      exp2<-min(round(log2(ny2))+1,amax);

      scale <- vector()

      j<-0;
      for (m in seq(exp1,exp2-1))
      {
        jj<-inter-1;
        for (n in seq(0,jj)) 
        {
          a<-2^(m+n/inter);
          j<-j+1;
          scale[j]<<-a;
        }
    }})
    } 
   

## now infers the corresponding periods
 ##  omega0<-1/2*(k0/aa+sqrt(2+k0*k0)/aa);
  omega0<-k0/scale;
  period<-1./omega0*2*pi*deltat;

  x<-y
  n<-length(x);

  k<-(0:(n-1))*2*pi/n

  f<-stats::fft(x);

  J<-length(scale);

  wave<-matrix(as.complex(0),nrow=n,ncol=J);

  if (calcmask) mask<-matrix(TRUE,nrow=n,ncol=J);
  if (deriv) dwave <- wave

  nn <-length(k);
  for(a1 in seq(along = scale)) {
     expnt<- -(scale[a1]*k - k0) ^2/2; ## psi_hat(a*k)
 ##   norm=sqrt(scale[a1]*k[2])*(pi^(-0.25))/sqrt(nn); 
    norm=2; 
    daughter=norm*exp(expnt);
    wave[,a1]=stats::fft(f*daughter ,inverse=TRUE)/n
    if (deriv) dwave[,a1]=stats::fft(f*daughter*(-(1i*k)) ,inverse=TRUE)/n;
    if (calcmask) {
      mask[1:min(n,ceiling(sqrt(2)*scale[a1])),a1] = NA;
      mask[(n-min(n,ceiling(sqrt(2)*scale[a1]))):n,a1] = NA;
    }
  }
    if (deriv)  {dwave <- dwave / (-1i*wave) }
    xx <-  as.array(xx);
    yy <- as.array(period)
    attr(wave,"time") <- xx
    if (calcmask) attr(wave,"mask") <- mask
    attr(wave,"period") <- yy
    attr(wave,"scale") <- scale
    attr(wave,"class") <- "wavelet"
    attr(wave,"wavelet") <- "morlet"
    attr(wave,"parameters") <- list(k0=k0,inter=inter,deltat=deltat)
    if (deriv) attr(wave,"deriv") <- dwave
    wave
}

#' @rdname cwt_morlet
#' @export
plot.wavelet <- function (x,resx=400,resy=300,xlab="Time",ylab="Period",scaling_correction=0,col=col_wavelet,legend=FALSE,Mode=Mod,plotMask=TRUE,...)

{
  require(fields)
  xx <- attr(x,"time")
  period <- attr(x,"period")
  mask <- attr(x,"mask")
  if (!plotMask) mask[,] = 1  # do not hide influence cone
  aa <- attr(x,"scale")
  thin_factor_xx <- max(ceiling(length(xx)/ resx),1)
  thin_factor_yy <- max(ceiling(length(aa)/ resy),1)
  subx <- seq(thin_factor_xx,length(xx),thin_factor_xx)
  suba <- seq(thin_factor_yy,length(aa),thin_factor_yy)
  wave_scaled <- Mode(x[subx,suba])*mask[subx,suba]
  for (i in 1:length(suba)) wave_scaled[,i] <- wave_scaled[,i]/(aa[i]^scaling_correction)
  if (legend) par(oma=c(2,2,2,5))
  image(xx[subx],period[suba],wave_scaled,log="y",ylab=ylab,xlab=xlab,axes=FALSE,col=col,...)
  axis(1)
  axis.log10(2,"")
  if (legend) {par(oma=c(2,2,2,2))
  image.plot(xx[subx],period[suba],wave_scaled,legend.only=TRUE,log="y",ylab=ylab,xlab=xlab,axes=FALSE,col=col,...)
  par(oma=c(2,2,2,5))}
}

#' @rdname cwt_morlet
#' @export
powerspectrum.wavelet <- function (wave,...)
{
  aa <- attr(wave,"scale")
  yy <- attr(wave,"period")
  xx <- attr(wave,"time")
  f <- 1./yy
  P <- vector("double",length(f))
for (j in seq_along(f)) P[j] <- mean(Mod(wave[,j])^2,na.rm=TRUE)
  data.frame("frequency"=f,"power"=P)
}

analytic.ridge <- function (wave,ridge_only = FALSE,plim=c(-Inf,+Inf))
{
  require("signal")
  p <- attr(wave,"period")
  P <- length(p)
  pout <- which(p<plim[1] | p>plim[2])
  s <- attr(wave,"scale")
  t <- attr(wave,"time")
  deltat <- diff(attr(wave,"time"))[1]
  WR <- wave
  WR[,pout] <- NA
  # scale for scalogram
  for (i in seq(1,P)) WR[,i] <- WR[,i]/sqrt(s[i])

  WP <- Arg(WR)
  for (i in seq(1,P) ) WP[,i] <- unwrap(WP[,i])
  for (i in seq(1,P) ) WP[,i] <- c(diff(WP[,i]),NA)/(2*pi)/deltat*p[i]
  WR[which(abs(WP-1) > 0.025)] <- NA
  if (ridge_only) 
  { 
   R <- vector()
   for (i in seq(1,length(t))) 
    {M <- which.max(Mod(WR[i,])) ; if (any(M)) R[i] <- WR[i,M] }
   R
  } else 
  WR

}
