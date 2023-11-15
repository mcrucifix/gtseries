
# compute correlation function for sin / cos basis
# generally small matrices need to be computed, at most 10x10
# so this is not the most critical area for optimisation

#' ancillary functions and definitions
pi2 <- pi*pi;
sinxx <- function(x) {if (x==0) return (1) else return (sin(x)/x) } 
fac <- function(x)  {return (sinxx(x) * pi2 / (pi2 - x*x))} ;
# fac <- function(x)  {return (cos(x))};
cossin <- function(x){if (x[2]==1) return(cos(x[1]*seq(0,(x[3]-1)))) else return(sin(x[1]*seq(0,(x[3]-1))))}

hann <- function(N) {return (1-cos(2*pi*seq(0,(N-1))/(N-1)))};


#' corr between two frequencies
#'
#'given two unitless frequencies (expressed as freq*(1/(N-1)), compute crossproduct of the two series
#'coscos, cossin, sincos, sinsin
#'as follows
#' / coscos   cossin \
#' |                 |
#' \ sincos   sinsin /
quarto <- function(f1,f2){
  Fplus <- fac(f1+f2);
  Fminus <- fac(f1-f2);
  cosminus <- 0.5*cos(f1-f2) * Fminus ;
  cosplus <- 0.5*cos(f1+f2) * Fplus;
  sinminus <- 0.5*sin(f1-f2) * Fminus;
  sinplus <- 0.5*sin(f1+f2) * Fplus;
  return ( matrix(c( cosplus+cosminus, 
                    sinplus+sinminus,
                    sinplus-sinminus,
                    cosminus-cosplus), 2,2) )
}

#' alpha matrix corresponding to a vector of frequencies
alpha <- function(Fs,N){
  Fn <- Fs*((N-1)/2);
  n <- length(Fs)
  DN <- 2*n;
  out <- matrix(0,DN,DN)
  for (i in seq(n)){
    for (j in seq(i,n)){
       Q <- quarto(Fn[i], Fn[j]);
       out[(2*i-1):(2*i) ,  (2*j-1):(2*j) ] <- Q
       out[(2*j-1):(2*j) ,  (2*i-1):(2*i) ] <- t(Q)
    }
  }
  if (Fs[1]==0) out=out[-2,-2]
  return(out)
}

#' construct a basis of frequencies Fs (expressed as omega 
basis <- function(Fs,N){
  n  <- length(Fs)
 return ( 
         apply( 
               cbind( as.vector(t(matrix(rep(Fs,2),n,2)))  , rep(c(0,1), n), rep(N,2*n))
               , 1, cossin) )
}

#' scalar product of two vectors weighted by HAnn window
sprod <- function(x,y){
  N <- length(x)
  return (  sum ( x * hann(N)*y ) / N  )
#   return (  sum ( x * y ) / N  )
}

#' orthonormalised basis based on existing basis and 
#' corresponding alpha matrix
gm_orth <- function(A) { 
  S <- svd(A) 
  sweep(S$v, 2, sqrt(S$d), '/')
}

calc_amp_fast <- function(x, nfreq, Fs){
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

     A <- alpha(c(0,Fs[1:j]),N)
     Coefs <- gm_orth(A)
     B <- xx %*% Coefs

#      B <- pracma::gramSchmidt(xx)$Q
#      Coefs = t(xx) %*% B


#      aa <- as.numeric(residu) %*% B
     aa <- t(apply(B,2, function(b) sprod(residu,b)))
#      aa <- sprod(residu, B)
     signal <- B %*% t(aa)
     residu <- residu - signal
  }
  
#   aa <- as.numeric(x) %*% B
     aa <- t(apply(B,2, function(b) sprod(x,b)))
  sol = 2.*solve(t(Coefs), t(aa))
  sol1 <- matrix(sol[-1],2, nfreq)
  
  Freqs <- c(0,Fs)
  amps <- c(sol[1], sqrt ( apply(sol1^2, 2, sum)  ))
  Phases <- -c(0, ( apply(sol1, 2, function(i) Arg(i[1] + 1i*i[2])  )))
  
  OUT = data.frame(freqs=Freqs, amps=amps, Phases=Phases) 
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
#' @importFrom pracma GramSchmidt
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
#' @export mfft_real_fast
mfft_real_fast <- function(xdata, nfreq=5){
  xdata = stats::as.ts(xdata)
  dt = deltat(xdata)
  startx = stats::start(xdata)[1]
  N = length(xdata)
  residu = xdata

  # will withold the definitive frequencies
  Fs <- rep(-1, nfreq)
    k=1
    OUT <- calc_amp_fast(xdata, nfreq, Fs)

    Fs <- OUT$freqs[-1]
    Freqs <- OUT$freqs
    amps <- OUT$amps
    Phases <- OUT$Phases

    print ('before correction')
    print (OUT)

    # correction 

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
    for (i in seq(nfreq+1)) xdata_synthetic = xdata_synthetic + OUT$amps[i]*cos(OUT$freqs[i]*seq(N) + OUT$Phases[i])
    OUT2 <- calc_amp_fast(xdata_synthetic, nfreq, Fs)
  
    print(OUT2)
    print(OUT$freqs - OUT2$freqs)
    # adjust to units
    OUT$freqs = OUT$freqs + (OUT$freqs - OUT2$freqs)
    OUT$amps = OUT$amps + (OUT$amps - OUT2$amps)
    OUT$Phases = OUT$Phases + (OUT$Phases - OUT2$Phases)
    OUT$freqs <- OUT$freqs/dt
    OUT$Phases <- OUT$Phases - startx*OUT$freqs
    return(OUT)
}

