#' Theoretical Spectrum of Autoregressive Process
#'
#' Calculates the theoretical density frequency of an
#' autoregressive process, according to equation
#' (3.3) of Akaike (1969).
#'
#' @param f numeric; vector of frequencies at which to calculate the spectrum.
#' @param ar numeric; vector of autoregressive coefficients.
#' @return Numeric vector or matrix of theoretical spectral density values.
#' @export
arspec <- function (f, ar) {
  Exp <- exp(-1i*2*pi*outer(f, seq_along(ar), "*"))
  1./abs(1-Exp%*%ar)^2
}
