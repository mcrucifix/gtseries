t <- seq(1024)

x_orig <- cos(t*0.54423167+0.00) + 1.3 * cos(t*0.441332+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) 
# x_orig <- 1.4*cos(t*0.54423167+0.41)
xdata <- x_orig

analyse <- function(xdata, nfreq){
  
  N <- length(xdata)
  hann <- function(N) {return (1-cos(2*pi*seq(0,(N-1))/(N-1)))};
  hN <- hann(N)
  t <- seq(N)
  power <- function(x) (Mod(fft(x))^2)[1:floor((N-1)/2)]
  
  hprod <- function(x,y) sum(x * hN * y)/N
  
  quarto <- function(f1,f2){
    M <- matrix(c(sum(cos(f1*t)*cos(f2*t)*hN), 
                 sum(cos(f1*t)*sin(f2*t)*hN), 
                 sum(sin(f1*t)*cos(f2*t)*hN), 
                 sum(sin(f1*t)*sin(f2*t)*hN)),2,2)
    M <- M/N
    return(M)
  }
  
  
  
    # B_i = sum[m<= i](A_im) f_m the orthogonal base signals
    # x_m what remains of the signal at time step m
  
  nu <- rep(0,nfreq)
  S <- rep(0,nfreq)
  nu <- rep(0,nfreq)
  phi <- rep(0,nfreq)
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
    fbase <- freqs[which.max(power(hx))]
    brackets <- c(fbase-pi/N, fbase+pi/N);
    thx <- t(hx)
    
    
    fmax = cmna::goldsectmax(function(f) 
                             (thx %*% cos(f*t))^2 + (thx %*% sin(f*t))^2, 
                             brackets[1], brackets[2], tol=1.e-10, m=9999)
    
    nu[m] = fmax
    
    # determine amplitude and phase
    
    # Be U = the matrix (ccoscos/cossin/sincos/sinsin) and its corresponbind crossproduct V
    # Be C the cos vector, and S the sin vector. They are _not_ orthogonal
    # Be X the signal for which we look for amplitude and phase
    # We are looking at the projection of U in the space spanned by C and S. This is
    
    q <- quarto(fmax, fmax)
    
    if (fmax > 1.e-10){
    xx <- rbind(cos(fmax*t), sin(fmax*t))
    prod <- xx %*% hx/N
    a  <- solve(q, prod)
    phase[m] <- -atan(a[2]/a[1])    
    } else {
      phase[m] = 0.
    }
    
    # phi[m] is < 
    
    f[[m]] <- cos(fmax*t + phase[m])
    phi[m] <- prod[1]*cos(phase[m]) - prod[2]*sin(phase[m])   # must be a way to be faster
    # should be equivalent to
    ## phi_test <- sum(f[[m]] * hx)/N    # to be verifed - -> ok
    
    
    for (i in seq(m)) Q[m,i] = hprod(f[[m]],f[[i]])  # so remember the convetion
                                                 # these are symmetric matrices
                                                 # and we only populate the j <= m 
    
    # before normalisation
    # B[[m]] = sum_j=1,m A[m,j]*f[j] et
    # B[[m]] = ( f[m] - sum_(1,m-1) (f[m] %*% B[i]) * B[i]
    # B[[m]] = ( f[m] - sum_(1,m-1) 
    #                (f[m] * sum(j=1,i A_j,i  f[i]) ) *
    #                               ( sum (A_j=1, i) A_j,i f[i] )
    
    #/ ||B[[m]]||
    
    # on peut alors verifier que
    
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
    
    # donc la on a, a ce stade, B[[m]] = sum [A[m,j]] * f[j]
    
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
    
    # j=1
    # amp
  #   [j]=0
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
  
  return(data.frame(amps=amp, freqs=nu, Phases=phase))
}
  
     
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
mfft_real_ter <- function(xdata, nfreq){
    N <- length(xdata)

    xdata = stats::as.ts(xdata)
    dt = deltat(xdata)
    startx = stats::start(xdata)[1]
#     print('startx')
#     print(startx)
    N <- length(xdata)
    OUT <- analyse(xdata, nfreq)

#     print('OUT before corr')
#     print(OUT)
#     Freqs <- OUT$freqs
#     amps <- OUT$amps
#     Phases <- OUT$Phases
    # correction 

    # correction  (methode 2)
#   }
   xdata_synthetic <- rep(0,N)
    for (i in seq(nfreq)) xdata_synthetic = xdata_synthetic + OUT$amps[i]*cos(OUT$freqs[i]*seq(N) + OUT$Phases[i])
    OUT2 <- analyse(xdata_synthetic, nfreq)
#     print ('Synthetic')
#     print(OUT2)
#     print ('Corrections')
#     print(OUT$freqs - OUT2$freqs)
    # adjust to units
    OUT$freqs = OUT$freqs + (OUT$freqs - OUT2$freqs)
    OUT$amps = OUT$amps + (OUT$amps - OUT2$amps)
    OUT$Phases = OUT$Phases + (OUT$Phases - OUT2$Phases)
    OUT$freqs <- OUT$freqs/dt
    OUT$Phases <- OUT$Phases - startx*OUT$freqs
#     print ('After corrections')
#     print(OUT)
    return(OUT)
}

