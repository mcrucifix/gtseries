
Q <- function(wT) {sin(wT)/(wT)*(pi^2)/(pi^2-(wT)^2)}
Qprime <- function(y) {ifelse(y==0, 0, pi^2/(pi^2-y^2)/y*(cos(y)+(sin(y)/y)*(3*y^2-pi^2)/(pi^2-y^2)))}
Qsecond0 <- 2/pi^2 - 1./3. 

analyse <- function(xdata, nfreq, fast = TRUE, nu = NULL){
  

  # nu can be provided, in which case the frequencies
  # are considered to be given
  
  N <- length(xdata)
  T <- N / 2
  hann <- function(N) {return (1-cos(2*pi*seq(0,(N-1))/(N)))};
  hN <- hann(N)
  t <- seq(N)-1
  power <- function(x) (Mod(fft(x))^2)[1:floor((N-1)/2)]
  
  hprod <- function(x,y) sum(x * hN * y)/N
  
  if (fast){
  quarto <- function(f1,f2){
    fm <- (f1-f2)*T
    fp <- (f1+f2)*T
    Qm <- ifelse(fm == 0, 1 ,  Q(fm))
    Qp <- ifelse(fp == 0, 1 ,  Q(fp))
    cosm <- cos(fm)*Qm
    cosp <- cos(fp)*Qp
    sinm <- sin(fm)*Qm
    sinp <- sin(fp)*Qp
    M <- 0.5 * matrix(c( 
            cosm + cosp , 
            sinm + sinp , 
            - sinm + sinp , 
            cosm - cosp ), 2 , 2 )
    return(M)
  }} else 
  {
  quarto <- function(f1,f2){
    M <- matrix(c(sum(cos(f1*t)*cos(f2*t)*hN), 
                 sum(cos(f1*t)*sin(f2*t)*hN), 
                 sum(sin(f1*t)*cos(f2*t)*hN), 
                 sum(sin(f1*t)*sin(f2*t)*hN)),2,2)
    M <- M/N
    return(M)
  }
  }
  
  
  
    # B_i = sum[m<= i](A_im) f_m the orthogonal base signals
    # x_m what remains of the signal at time step m
  
  S <- rep(0,nfreq)
  if (is.null(nu)) nu <- rep(NA,nfreq)
  phase <- rep(0,nfreq)
  amp <- rep(0,nfreq)
  A <- matrix(0,nfreq,nfreq)
  Q <- matrix(0,nfreq,nfreq)
  f <-  list()
  x <-  list()
  B <- list()
  # S[m] = hprod (f[m], B[m])
  freqs = 2.*pi*seq(0,(N-1))/N
  x[[1]] = xdata
  
   
  for (m in seq(nfreq)){
  hx <- hN*x[[m]]
    if (is.na(nu[m])){
    # are frequencies already provided ? 
    fbase <- freqs[which.max(power(hx))]
    brackets <- c(fbase-pi/N, fbase+pi/N);
    thx <- t(hx)
    
    # after profiling, the fastest seems the first option below
    # this is the weak link
    # coding the function in c might be an option
    
    fmax = cmna::goldsectmax(function(f) 
                             {
                              ft <- f*t
                             (thx %*% cos(f*t))^2 + (thx %*% sin(f*t))^2}, 
                             brackets[1], brackets[2], tol=1.e-10, m=9999)
    
   #  fmax = cmna::goldsectmax(function(f) {
   #                           ft <- f*t
   #                           (sum(hx * cos(ft)))^2 + (sum(hx %*% sin(ft)))^2}, 
   #                           brackets[1], brackets[2], tol=1.e-10, m=9999)

   #  fmax = cmna::goldsectmax(function(f) {
   #                           ft <- f*t
   #                           Mod(sum(hx * exp(1i*ft)))}, 
   #                           brackets[1], brackets[2], tol=1.e-10, m=9999)



    nu[m] <- fmax
} else {
   fmax <- nu[m] }   # else, use the provided frequency
    
    # determine amplitude and phase
    
    # Be U = the matrix (ccoscos/cossin/sincos/sinsin) and its corresponbind crossproduct V
    # Be C the cos vector, and S the sin vector. They are _not_ orthogonal
    # Be X the signal for which we look for amplitude and phase
    # We are looking at the projection of U in the space spanned by C and S. This is
    
    q <- quarto(fmax, fmax)
    if (fmax > freqs[2]/2){
    xx <- rbind(cos(fmax*t), sin(fmax*t))
    prod <- xx %*% hx/N

    # to do : given that q is only 2x2 we do not need the solve function
    # it is pretty easy to o the diagonalisation by hand. 
    a  <- solve(q, prod)
    phase[m] <- -atan(a[2]/a[1])    
    } else {
      phase[m] = 0.
      nu[m] = 0.
    }
    
    f[[m]] <- cos(fmax*t + phase[m])
    
    if (fast)  # we use analytical forms of the products
    {
       
    for (i in seq(m))
    {

        num <- (nu[m] - nu[i])*T
        nup <- (nu[m] + nu[i])*T
        phim <- (phase[m] - phase[i])
        phip <- (phase[m] + phase[i])

        Qm <- ifelse(num == 0, 1 ,  Q(num))
        Qp <- ifelse(nup == 0, 1 ,  Q(nup))
        cosm <- cos(phim) * cos(num) * Qm
        sinm <- sin(phim) * sin(num) * Qm
        cosp <- cos(phip) * cos(nup) * Qp
        sinp <- sin(phip) * sin(nup) * Qp

        Q[m,i] = 0.5 * (cosm + cosp - sinm - sinp)
    }
    } else {
    for (i in seq(m)) 
    {Q[m,i] = hprod(f[[m]],f[[i]])  # so remember the convetion
                                                 # these are symmetric matrices
                                                 # and we only populate the j <= m 
    }
    }
    

    #
    # before normalisation
    # B[[m]] = sum_j=1,m A[m,j]*f[j] et
    # B[[m]] = ( f[m] - sum_(1,m-1) (f[m] %*% B[i]) * B[i]
    # B[[m]] = ( f[m] - sum_(1,m-1) 
    #                (f[m] * sum(j=1,i A_j,i  f[i]) ) *
    #                               ( sum (A_j=1, i) A_j,i f[i] )
    
    #/ ||B[[m]]||
    
    # one can then verify that:
    
    # B[[m]] %*% B[j]] = 
    # ( f[m] - sum_(1,m-1) (f[m] %*% B[i]) * B[i] ) %*%  B[j] =
    # ( f[m] - (f[m] %*% B[j]) * B[j] ) %*%  B[j] =
    # ( f[m] %*% B[j] - f[m] %*% B[j]) = 0 
    
    # before normalisation
    A[m,] = 0
    A[m, m] = 1.
    if (m>1){
      # fmbi = f[m] %*% B[i]
      fmbi <- rep(0,(m-1))
      for (j in seq(m-1)) for (i in seq(j)) fmbi[j] = fmbi[j] + A[j,i]*Q[m,i]
      for (j in seq(m-1)) for (i in seq(j,(m-1))) A[m,j] = A[m,j] - fmbi[i]*A[i,j]
    }
    
    # so, at this stage, we have B[[m]] = sum [A[m,j]] * f[j]
    
    # B[[2]] = (A[2,1]*f[1] + A[2,2]*f[2])
    # B[[3]] = (A[3,1]*f[1] + A[3,2]*f[2] + A[3,3] * f[3] ) 
    
    norm = 0
    if (m > 1){
    for (i in seq(m)) {
      norm = norm + (A[m,i]*A[m,i]*Q[i,i])
      if (i>1) for (j in seq(i-1)) norm  = norm + 2*A[m,i]*A[m,j]*Q[i,j]
    }
    } else {
      norm = A[m,m]*A[m,m]*Q[m,m]
    }
    
    
    A[m,] = A[m,] / sqrt(norm)
      
    # S[m] = x_orig %*% B[[m]] = sum hprod(x[[1]],f[[j]])*A[m,j]
     
    S[m] = 0. 
    for (j in seq(m))  S[m] = S[m] + hprod(x[[1]],f[[j]])*A[m,j]
    # not necessary, for verification only
    B[[m]]=0
    for (j in seq(m)) B[[m]] = B[[m]] + A[m,j] * f[[j]]
    
    # if you are curious, the amplitude of the signal will be
    # amp[j] = contribution of the sum( S[m]*B[m]) to the f[[m]]
    # amp[j] = sum(m in seq(j)) S(m) * A[m,j]
    
    # for (i in seq(m)) amp[j] = amp[j]+S[m]*A[m,j]
    # and the residual signal
    
    # f[[m+1]] = f[[m]] - S[m] * B[[m]]

    x[[m+1]] = x[[m]]
    for (j in seq(m))  x[[m+1]] = x[[m+1]] - S[m] * A[m,j]*f[[j]]
  
  }
  
  # amplitudes
  
  for (m in seq(nfreq)){
    amp[m]=0;
    for (j in seq(m)) amp[m] = amp[m] + A[m,j]*S[m]
  }
  
  OUT = data.frame(nu=nu, amp=amp, phase=phase) 
  return(OUT)
}
  
     
#' Modified Fourier transform  for real series (variant)
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
#' @param xdata The data provided either as a time series (advised), or as a vector. 
#' @param nfreq is the number of frequencies returned, must be smaller that the length of  xdata.
#' @param fast (default = TRUE) uses analytical formulations for the crossproducts involving sines and cosines
#' @param correction:  0: no frequency correction (equivalent to Laskar); 1 : frequency correction using linear approximation ; 2: frequency correction using sythetic data
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
#' @export mfft_real_ter
  # will withold the definitive frequencies
mfft_real_ter <- function(xdata, nfreq=5, correction=1, fast=TRUE){
    N <- length(xdata)
    N2 <- N/2.
    xdata = stats::as.ts(xdata)
    dt = deltat(xdata)
    startx = stats::start(xdata)[1]
    N <- length(xdata)
    OUT <- analyse(xdata, nfreq, fast)


    # correction  (methode 2)
#   }
  if (correction == 2){
   xdata_synthetic <- rep(0,N)
   t <- seq(N)-1
   for (i in seq(nfreq)) xdata_synthetic = xdata_synthetic + OUT$amp[i]*cos(OUT$nu[i]*t + OUT$phase[i])
   OUT2 <- analyse(xdata_synthetic, nfreq, fast)
   OUT$nu = OUT$nu + (OUT$nu - OUT2$nu)
   OUT$amp = OUT$amp + (OUT$amp - OUT2$amp)
   OUT$phase = OUT$phase + (OUT$phase - OUT2$phase)
   } else if (correction == 1){

    for (j in seq(nfreq)){
       epsilon = OUT$amp[j] * Qprime(-2 * OUT$nu[j] * N2)*cos(2 * OUT$nu[j] * N2 + 2 * OUT$phase[j])
       print ('epsilon')
       print (epsilon)
       if ((j+1) <= nfreq) { for (s in seq(j+1, nfreq)) {
        epsilon = epsilon + OUT$amp[s]  * 
         ( 
          Qprime( (OUT$nu[s] - OUT$nu[j])*N2)*cos((OUT$nu[j] - OUT$nu[s])*N2 + OUT$phase[j] - OUT$phase[s] ) -
          Qprime(( OUT$nu[s] + OUT$nu[j])*N2)*cos((OUT$nu[j] + OUT$nu[s])*N2 + OUT$phase[j] + OUT$phase[s] ) )
     }}
    epsilon = epsilon / Qsecond0 / N2 / OUT$amp[j]

    OUT$nu[j] = OUT$nu[j] - epsilon
    }
    OUT <- analyse(xdata, nfreq, fast, nu = OUT$nu)
  }

  # account for tseries scaling (i.e. uses the vaue of deltat and start encoded
  # in the time series object 
    OUT$nu <- OUT$nu / dt
    OUT$phase <- OUT$phase - startx*OUT$nu

  # rearrange terms to avoid negative amplitudes and have phases in the right quandrant
  # these are cosines so even functions 

    to_be_corrected <- which (OUT$amp < 0)
    if (length(to_be_corrected)){
      OUT$amp[to_be_corrected] <- - OUT$amp[to_be_corrected]
      OUT$phase[to_be_corrected] <- OUT$phase[to_be_corrected] + pi 
    }

    to_be_corrected <- which (OUT$nu < 0)
    if (length(to_be_corrected)){
      OUT$nu[to_be_corrected] <- - OUT$nu[to_be_corrected]
      OUT$phase[to_be_corrected] <- - OUT$phase[to_be_corrected] 
    }

   # finally, order termes by decreasing amplitude 
    o <- order(OUT$amp, decreasing = TRUE)
    OUT$amp <- OUT$amp[o]
    OUT$nu <- OUT$nu[o]
    OUT$phase <- OUT$phase[o]

    OUT$phase <- (OUT$phase + (2*pi)) %% (2*pi)


    # rename for class compatibility
    names(OUT) <- c("Freq","Amp","Phases")
    attr(OUT, "class") <- "mfft_deco"
    attr(OUT, "data")  <- xdata
    attr(OUT, "nfreq")  <- nfreq
    return(OUT)
}

