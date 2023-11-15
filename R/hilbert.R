#' Hilbert extension
#'
#' Computes the extension of a real valued signal to a complex signas with
#' the Hilbert transform
#'
#'
#' @param x Input array
#' @return Analytic signal, of length \code{n}, returned as a complex vector or
#'
#' @export
hilbert_extension <- function(x) {

  n = length(x)
  # construct multiplication vector
  if (n %% 2 == 0) {
    v <- c(0, rep(-1i, n / 2 - 1), 1, rep(1i, n / 2 - 1))
  } else {
    v <- c(0, rep(-1i, (n - 1) / 2), rep(1i, (n - 1) / 2))
  }

  # compute the Hilbert transform
  return( x + 1i* Re(stats::fft(stats::fft(x)*v, inverse=TRUE))/n )
}
