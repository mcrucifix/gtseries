
pisquare <- pi^2

Q <- function(wT) {sin(wT)/(wT)*(pisquare)/(pisquare-(wT*wT))}
Qprime <- function(y) {ifelse(y==0, 0, pisquare/(pisquare-y*y)/y*(cos(y)+(sin(y)/y)*(3*y*y-pisquare)/(pisquare-y*y)))}
Qsecond0 <- 2/pisquare - 1./3. 

mfft_analyse <- function(xdata, nfreq, fast = TRUE, nu = NULL, minfreq=NULL, maxfreq=NULL){
  
  # nu can be provided, in which case the frequencies
  # are considered to be given
  if (!is.null(nu))
  {
    # we will need two versions of nfreq, one for the cosine and one for the cosine
    # if if there is a zero frequency we only need it once
    nutemp <- unlist ( lapply( nu, function(a) if (a==0) return (a) else return (  c(a,-a)))  )
    phase <-  unlist ( lapply( nu, function(a) if (a==0) return (0) else return (  c(0, pi/2)))  )
    nu <- nutemp
    if (length(nu) < 2*nfreq) {
      # this will typically occur if one zero frequency was provided
      nu[2*nfreq] <- NA
      phase[2*nfreq] <- NA
    }
  } else 
  {
    nu <- rep(NA,2*nfreq)
    phase <- rep(NA,2*nfreq)
  }

  N <- length(xdata)
  T <- N / 2
  hann <- function(N) {return (1-cos(2*pi*seq(0,(N-1))/(N)))};
  hN <- hann(N)
  t <- seq(N)-1
  power <- function(x) (Mod(fft(x))^2)[1:floor((N-1)/2)]
  
  hprod <- function(x,y) sum(x * hN * y)/N
  
  # B_i = sum[m<= i](A_im) f_m the orthogonal base signals
  # x_m what remains of the signal at time step m
  
  S <- rep(0,2*nfreq)
  Prod <- rep(0,2*nfreq)
  amp <- rep(NA,2*nfreq)
  A <- matrix(0,2*nfreq,2*nfreq)
  Q <- matrix(0,2*nfreq,2*nfreq)
  f <-  list()
  x <-  list()
  B <- list()
  # S[m] = hprod (f[m], B[m])
  freqs = 2.*pi*seq(0,(N-1))/N
  x[[1]] = xdata
  
   
  for (m in seq(2*nfreq)){
  hx <- hN*x[[m]]
    # there is a little tweak here: 
    # if we have reached 2*nfreq and not nu[m] is na, then
    # it means that we have encountered a constant somewhere, and that last step should be skipped. 
    if ( ! ( m == 2*nfreq && is.na(nu[m]))) {
    # ok we are on business
    if (is.na(nu[m])){
    # is the "m" frequency already provided ? 
    # either there were provided from the start and developped  in the first lines of this routine
    # or either they were not provided, but we are in this case where it was set in "m-1" because
    # we identified a non-null frequency
    # in both configurations, we look at a frequency with the fourier transform
     fbase <- freqs[which.max(power(hx))]
     brackets <- c(fbase-pi/N, fbase+pi/N);
    # and then further corrects the brackets if minfreq and maxfreq where provided
     brackets[1] <- max(minfreq, brackets[1])
     brackets[2] <- min(maxfreq, brackets[2])
     thx <- t(hx)
    
    # after profiling, the fastest seems the first option below
    # this is the weak link
    # coding the function in c might be an option
    
    #
     tomax <-  function(t) {
       function(f) {
       ft <- f*t
       a <- hx %*% cbind(cos(ft), sin(ft))
       a[1]*a[1] + a[2]*a[2]
     } 
     }


#     print ("x[[m]][1,2]")
#     print (as.double(x[[m]][1:2]))

    localxdata <- as.double(x[[m]])
#     print ('call it')
#    OUT <- .C("fastgsec", as.double(minfreq), as.double(maxfreq), 
#              as.integer(N), localxdata , as.double(rep(0,N)), 
#              outfreq = 0., DUP=TRUE)
#     print ('called it')

#     plot (Mod(fft(x[[m]])[0:(N/2)])^2 , type='l')

   fmax = cmna::goldsectmax(tomax(t), 
                            brackets[1], brackets[2],
                            tol=1.e-10, m=9999)
   

    # fmax = cmna::goldsectmax(function(f) {
    #                          ft <- f*t
    #                          (sum(hx * cos(ft)))^2 + (sum(hx %*% sin(ft)))^2}, 
    #                          brackets[1], brackets[2], tol=1.e-10, m=9999)

   #  fmax = cmna::goldsectmax(function(f) {
   #                           ft <- f*t
   #                           Mod(sum(hx * exp(1i*ft)))}, 
   #                           brackets[1], brackets[2], tol=1.e-10, m=9999)


    #print (sprintf("fmaxlocal: %.2f ; fmax with C: %.2f", 
    #               fmax, OUT$outfreq))

#    fmax <- OUT$outfreq

    if (fmax > freqs[2]/2){
      # if we really identified a frequency
      # we are in this case where the frequency was not a priori provided
      # we se set this phase to zero, and the next one to pi/2, and also set the frequencies accordingly
    phase[m] <- 0.
    nu[m]  <- fmax
    phase[m+1] <- pi/2. 
    nu[m+1] <- -fmax
    } else {
      # again we are still in thes case where a freqency was not a priori provided
      # but the more particular case where the frequency is (with tolerance) considered as zero
      # f[[m]] will effectively be constant. 
      phase[m] <- 0. 
      nu[m] <- 0.
    }}
    
    f[[m]] <- cos(nu[m]*t + phase[m])

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
      if (i>1) for (j in seq(i-1)) { norm  = norm + 2*A[m,i]*A[m,j]*Q[i,j]
    }}
    } else {
      norm <- A[m,m]*A[m,m]*Q[m,m]
    }
 
 
    
    A[m,] = A[m,] / sqrt(norm)
      
     
    Prod[m] = hprod(x[[1]],f[[m]])
    # les prods sont corrects

    S[m] = 0. 
    for (j in seq(m))  S[m] = S[m] + Prod[j]*A[m,j]

    # for (j in seq(m))  S[m] = S[m] + A[m,j] * hprod(x[[1]], f[[j]])

    # les Sm sont les projections dans la nouvelle base
    # not necessary, for verification only
    # computes the B and verify they are orthonormal
    # B[[m]]=0
    # for (j in seq(m)) B[[m]] = B[[m]] + A[m,j] * f[[j]]


    #print ('S[m] computed two different ways')
    #print (m)
    #print (S[m])
    #print (hprod(x[[1]], B[[m]]))
    #print ('end test')

    # for (j in seq(m)) print (sprintf("B%i%i : %3g", i, j, hprod(B[[j]], B[[m]])))
    # if you are curious, the amplitude of the signal will be
    # amp[j] = contribution of the sum( S[m]*B[m]) to the f[[m]]
    # amp[j] = sum(m in seq(j)) S(m) * A[m,j]
    
    # for (i in seq(m)) amp[j] = amp[j]+S[m]*A[m,j]
    # and the residual signal
    
    # x[[m+1]] = x[[m]] - S[m] * B[[m]]


    x[[m+1]] = x[[m]]
    for (j in seq(m))  x[[m+1]] = x[[m+1]] - S[m] * A[m,j]*f[[j]]
    #x[[m+1]] = x[[m+1]] - S[m] * B[[m]]
#     print('residu final')
#     print(m)
#     print (x[[m+1]][seq(20)])
  
  }

  }
  
  #xtest <- 0 
  #xtest2 <- 0 
  #coeftests <- rep(0,m)
  #phi1 <- 0.4 ; phi2 <- 0.1234
  # coefreals <- c(cos(phi1), -sin(phi1), cos(phi2), -sin(phi2))
  # for (j in seq(m)) xtest <- xtest + coefreals[j] * f[[j]]
  # for (j in seq(m)) for (i in seq(j)) xtest2 <- xtest2 + S[j] * A[j,i] * f[[i]]
  # for (j in seq(m)) for (i in seq(j)) coeftests[i] <- coeftests[i] + S[j] * A[j,i]
  # print ((xtest - x[[1]])[seq(20)])
  # print ((xtest2 - x[[1]])[seq(20)])
  # print ("cofreals")
  # print (coefreals)
  # print ("coftests")
  # print (coeftests)
  
  # donc: j'ai demontre que les deux xtests sont les memes
  # coefreals est la bonne solution
  # et je sais ussi que B = sum A_mj f_j
  # donc je devroais simplement avour les Sj A_mj
  

  # at htis point I would like to make a littel check
  # make A full
  # is A Q t(A) = BtB orthonormal ? 
  # for (m in seq(2,2*nfreq)) for (j in seq(m-1)) A[j,m]=A[m,j]
  # for (m in seq(2,2*nfreq)) for (j in seq(m-1)) Q[j,m]=Q[m,j]
  # print (A*t(Q)*t(A))
  # print ('END TEST')

  #print (A) 
  #print (A %*% S)
  # amplitudes
  
  mmax <- 2*nfreq
  if (is.na(nu[mmax])) mmax <- mmax - 1

  amp[1:mmax] = 0;
  for (m in seq(mmax)) for (j in seq(m)) amp[j] <- amp[j] + S[m] * A[m,j]

  amp2m <- amp;

  for (m in seq(mmax)){
    # if the perivous "m" was already of the same frequency, we can merge amplitudes and set phase
    if ((m > 1) && (nu[m-1] == -nu[m])){
        
        phase[m] <- Arg (-1i * amp[m] + amp[m-1])
        amp[m] <- sqrt(amp[m-1]^2 + amp[m]^2)
        amp[m-1] <- NA
        phase[m-1] <- NA
    }
  }
  
  # at this point we should have exactly nfreq non na values
  # we print a message if this is not right, and in that case I suspect some debugging will be needed. 

  valid_frequencies <- which(!is.na(amp))
  nu <- -nu[valid_frequencies]
  amp <- amp[valid_frequencies]
  phase <- phase[valid_frequencies]
  if (length(valid_frequencies) != nfreq) message (sprintf("something goes wrong : %i valid frequencies, and nfreq = %i", valid_frequencies, nfreq))

  OUT = data.frame(nu=nu, amp=amp, phase=phase) 
  return(OUT)
}
  
     
#' Frequency Modified Fourier transform  for real series 
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
#' @param minfreq,maxfreq If provided, bracket the frequencies to be probed. Note this are
#'        more precisely angular velocities (2\pi / period), expressed in time-inverse units
#'        with the time resolution encoded in `xdata` if the latter is a time series. 
#' @param correction:  0: no frequency correction (equivalent to Laskar); 1 : frequency correction using linear approximation ; 2: frequency correction using sythetic data (both are documented in the Sidlichovsky and Nesvorny paper. 
#' @param nfreq is the number of frequencies returned, must be smaller that the length of  xdata.
#' @param fast (default = TRUE) uses analytical formulations for the crossproducts involving sines and cosines. 
#'        note: this is not really faster because the bottleneck is actually the goden section search. But more elegant. 
#' @return a `mfft_deco` object, based on a data.frame with columns "Freq", "Ampl" and "Phases". 
#' @author Michel Crucifix
#' @references
#' \insertRef{sidlichovsky97aa}{gtseries}
#' @examples
#' 
#' data(harmonic_sample)
#' spectrum <- mfft_real(harmonic_sample$data)
#' print(spectrum)
#'
#' @export mfft_real
mfft_real <- function(xdata, nfreq=5,  minfreq=NULL, maxfreq=NULL, correction = 1 , fast=TRUE){

  N <- length(xdata)
  N2 <- N/2.
  xdata = stats::as.ts(xdata)
  dt = deltat(xdata)

  my_minfreq <- ifelse (is.null(minfreq), 0, minfreq*dt)
  my_maxfreq <- ifelse (is.null(maxfreq), pi, maxfreq*dt)

  startx = stats::start(xdata)[1]
  N <- length(xdata)
  OUT <- mfft_analyse(xdata, nfreq, fast, NULL, my_minfreq, my_maxfreq)


  # correction  (methode 2)
  if (correction == 2){
    xdata_synthetic <- rep(0,N)
    t <- seq(N)-1
    for (i in seq(nfreq)) xdata_synthetic = xdata_synthetic + OUT$amp[i]*cos(OUT$nu[i]*t + OUT$phase[i])
    OUT2 <- mfft_analyse(xdata_synthetic, nfreq, fast, NULL, my_minfreq, my_maxfreq)
    OUT$nu = OUT$nu + (OUT$nu - OUT2$nu)
    OUT$amp = OUT$amp + (OUT$amp - OUT2$amp)
    OUT$phase = OUT$phase + (OUT$phase - OUT2$phase)
  } else if (correction == 1){

    for (j in seq(nfreq)){
      epsilon = OUT$amp[j] * Qprime(-2 * OUT$nu[j] * N2)*cos(2 * OUT$nu[j] * N2 + 2 * OUT$phase[j])
      if ((j+1) <= nfreq) {
        for (s in seq(j+1, nfreq)) {
          epsilon = epsilon + OUT$amp[s]  * 
           ( 
            Qprime( (OUT$nu[s] - OUT$nu[j])*N2)*cos((OUT$nu[j] - OUT$nu[s])*N2 + OUT$phase[j] - OUT$phase[s] ) -
            Qprime(( OUT$nu[s] + OUT$nu[j])*N2)*cos((OUT$nu[j] + OUT$nu[s])*N2 + OUT$phase[j] + OUT$phase[s] ) 
           )
        }
      }
      epsilon = epsilon / Qsecond0 / N2 / OUT$amp[j]

      OUT$nu[j] = OUT$nu[j] - epsilon
    }
    OUT <- mfft_analyse(xdata, nfreq, fast, nu = OUT$nu, my_minfreq, my_maxfreq)
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

