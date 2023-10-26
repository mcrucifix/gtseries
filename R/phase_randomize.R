#' Phase Randomizer
#'
#' randomizes the phases of a real signal, for non-linear time series analysis
#' It works by computing the FFT of the signal, then randomize the phases of all
#' complex numbers, but preserving the complex conjugate relationship guarantee
#' that the inverse transform will be real
#' @param  x the input time series
#' @importFrom stats runif
#' @return x the output time series. 
#'
phase_randomize <- function (x)
{
N <- length(x)

X <- fft(x)

k <- seq(2,floor(N/2)+1)

X[k] <- Mod(X[k])*exp(2*pi*1i*stats::runif(length(k)))
X[N-k+2] <- Conj(X[k])

if (N == 2*floor(N/2)) {k <- N/2+1
                        X[k] <- Mod(X[k])}

y <- fft(X,inverse=TRUE)
Re(y)/N
}

