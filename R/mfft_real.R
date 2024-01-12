calc_amp <- function(x, nfreq, Fs){
  N <- length(x)
  t <- seq(N)
  xx = matrix(rep(1,N),N,1)
  freqs = 2.*pi*seq(0,(N-1))/N
  power <- function(x) (Mod(fft(x))^2)[1:floor((N-1)/2)]
  residu <- x  - mean(x)
  for (j in seq(nfreq)) {
 
# temps de calcul peut etre gagne en C


    if (Fs[j] == (-1) ){
       hresidu = gsignal::hann(N)*residu
       fbase <- freqs[which.max(power(hresidu))]
       brackets <- c(fbase-pi/N, fbase+pi/N);
       thresidu <- t(hresidu)
       fmax = cmna::goldsectmax(function(f) 
                          (thresidu %*% cos(f*t))^2 + (thresidu %*% sin(f*t))^2, 
                          brackets[1], brackets[2], tol=1.e-10, m=9999)
       Fs[j] = fmax
    } 
 
     xx <- cbind(xx, cos(Fs[j]*t), sin(Fs[j]*t))
     
# cette partie doit pouvoir etre acceleree en tirant partit du fait
# que l'expression analytitque de xx%*xx est connue
# voir non-working code fmft_real.C
# Il faut utiliser l'integrale connue int_0^2pi cos(x)hann(x) 
# on utilise: cos(a)cos(b) = 0.5*(cos(a+b) - cos(a-b)) etc pour toutes les comb. cossin
# et int_0^2T cos(wx) = 1/W(sin[2T]-sin[0]) , ou si on applique le produit scalaire avec fenetre Hann
# et int_0^2T cos(wx)Hann[x] = (pi^2)(pi^2-w^2T^2)(sin(wT)/wT) * (coswT)
# et int_0^2T sin(wx)Hann[x] = (pi^2)(pi^2-w^2T^2)(sin(wT)/wT) * (sin(wT))
 # sin(wT)/(wT) = 1
 # et on se rappelle que la decomposition gramschmidt
 # c'est simplement, avec x_i le ie vecteur
 # xn_1 = x_1 / ||x_1||  
 # xn_(i+1) = x_(i+1) - sum_(j=1,i) x_(i+1)%*%x(i)/||x_(i+1)||
# en C, cela doit reduire le temps de calcul considerablement


     B <- pracma::gramSchmidt(xx)$Q
     Coefs = t(xx) %*% B
     aa <- as.numeric(residu) %*% B

     signal <- B %*% t(aa)
     residu <- residu - signal
  }
  
  aa <- as.numeric(x) %*% B
  sol = solve(t(Coefs), t(aa))
  sol1 <- matrix(sol[-1],2, nfreq)
  
  Freqs <- c(0,Fs)
  amps <- c(sol[1], sqrt ( apply(sol1^2, 2, sum)  ))
  Phases <- -c(0, ( apply(sol1, 2, function(i) Arg(i[1] + 1i*i[2])  )))
 

  OUT = data.frame(Freq=Freqs, Amp=amps, Phase=Phases) 
  attr(OUT,"class") = "mfft_deco"
  return(OUT)
}



#' Modified Fourier transform  for real series
#'
#' R-coded version of the Modified Fourier Transform
#' with frequency  correction, adapted to R. 
#' much slower than mfft (for complex numbers) as the latter is
#' mainly written in C, but is physically
#' more interpretable if signal is real, because
#' it is designed to have no imaginary part in the residual
#' A C-version should be supplied one day. 
#'
#' @importFrom cmna goldsectmax
#' @importFrom pracma gramSchmidt
#' @importFrom gsignal hann
#' @param xdata The data provided either as a time series (advised), or as a vector. 
#' @param nfreq is the number of frequencies returned, must be smaller that the length of  xdata.
#' @author Michel Crucifix
#' @references
#' \insertRef{sidlichovsky97aa}{gtseries}
#' @examples
#'
#' set.seed(12413)
#' t = seq(1024)
#' x_orig = cos(t*0.13423167+0.00) + 1.3 * cos(t*0.119432+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) + rnorm(1024)*0.12
#' OUT <- mfft_real(x_orig)
#' print(OUT)
#'
#' @export mfft_real
mfft_real <- function(xdata, nfreq=5){
  xdata = stats::as.ts(xdata)
  dt = deltat(xdata)
  startx = stats::start(xdata)[1]
  N = length(xdata)
  times <- (seq(N)-1)*dt+startx
  residu = xdata

  # will withold the definitive frequencies
  Fs <- rep(-1, nfreq)
    k=1
    OUT <- calc_amp(xdata, nfreq, Fs)

    Fs <- OUT$Freq[-1]


##   correction  METHOD ANALYTICAL. CURRENTLY NOT USED BUT SHOULD STILL BE EXPLORED
##  -------------------------------------------------------------------------------- 
##    Freqs <- OUT$Freq
##    amps <- OUT$Amp
##    Phases <- OUT$Phase
##     if (k < (nfreq)){
##       tau = (N-1) / 2.
##       Qp <- function(y) {
##         (1 / y) * (pi^2/ (pi^2-y^2))*(cos(y) + sin(y) / y * (3. * y^2 - pi^2)/(pi^2 - y^2))}
##       Q0 <- {2. / (pi^2) -1./3.}
##       
##       epsilon = rep(0,nfreq)
##       
## #       tau = (N-1) / 2.
##      
##       for (j in seq(k,nfreq-1)){
##          jj = j+1
##          ysj =(Freqs[(jj+1):(nfreq+1)]-Freqs[jj])*tau
##          epsilon[j] = sum(amps[(jj + 1):(nfreq+1)]*Qp(ysj)/
##                           (amps[jj] * Q0 * tau) * cos(ysj + Phases[(jj+1):(nfreq+1)]-Phases[jj]))
##       }
##       Epsilon=c(0,epsilon)
##       Freqs <- Freqs - Epsilon
##       Fs <- Freqs[2:(nfreq+1)]
##     }
## 

   
    # correction  (methode 2)

#   }
   xdata_synthetic <- rep(0,N)
     Fs <- rep(-1, nfreq)
    for (i in seq(nfreq+1)) xdata_synthetic = xdata_synthetic + OUT$Amp[i]*cos(OUT$Freq[i]*seq(N) + OUT$Phase[i])
    OUT2 <- calc_amp(xdata_synthetic, nfreq, Fs)
  
#     print(OUT2)
#     print(OUT$Freq - OUT2$Freq)
    # adjust to units
    OUT$Freq = OUT$Freq + (OUT$Freq - OUT2$Freq)
    OUT$Amp = OUT$Amp + (OUT$Amp - OUT2$Amp)
    OUT$Phase = OUT$Phase + (OUT$Phase - OUT2$Phase)
    OUT$Freq <- OUT$Freq/dt
    OUT$Phase <- OUT$Phase - startx*OUT$Freq

    O <- order(OUT$Amp, decreasing=TRUE)

    OUT$Amp = OUT$Amp[O]
    OUT$Freq = OUT$Freq[O]
    OUT$Phase = OUT$Phase[O]

  

    attr(OUT,"class") = "mfft_deco"
    attr(OUT,"nfreq") = nfreq
    attr(OUT,"xdata") = xdata

    total_ssq <- sum(xdata^2)

    return(OUT)
}

#' MFFT reconstruction
#' @rdname mfft_deco
#' @export reconstruct_mfft
reconstruct_mfft  <- function(M){
 if (!(attr(M,"class")=="mfft_deco")) stop ("object is not a MFFT decomposition")
 xdata <- attr(M,"xdata")
 nfreq <- attr(M,"nfreq")
 times <- seq(length(xdata))*dt(xdata) + startx(xdata)
 reconstructed <- lapply(seq(nfreq), function(i) M$Amp[i]*cos(M$Freq[i]*times + M$Phase[i]) )
}

#' MFFT ANOVA
#' @rdname mfft_deco
#' @export mfft_anova
reconstruct_mfft  <- function(M){
 if (!(attr(M,"class")=="mfft_deco")) stop ("object is not a MFFT decomposition")
 xdata <- attr(M,"xdata")
 nfreq <- attr(M,"nfreq")
 N <- length(xdata)
 times <- seq(length(xdata))*deltat(xdata) + start(xdata)
 reconstructed <- sapply(seq(nfreq), function(i) M$Amp[i]*cos(M$Freq[i]*times + M$Phase[i]) )
 cum_reconstruct <- apply(reconstructed, 1, cumsum)
 residual_vars <- apply(apply(cum_reconstruct, 1, function(i) xdata-i) , 2, function(j) sum(j^2))
 var0 <- sum(xdata^2)
 master_vars <- c(var0, residual_vars[-length(residual_vars)])
 p2s <- 2*seq(nfreq)
 p1s <- c(0, p2s[-length(p2s)])
 F <- (master_vars - residual_vars)/residual_vars * (N - p2s)/(p2s-p1s)
}





#' @rdname mfft_deco
#' @export
plot.mfft_deco <- function (M,...){
  O <- order(M$Freq)
  plot(M$Amp[O], M$Freq[O],'h',...)
}



