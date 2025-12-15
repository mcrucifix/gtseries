#' Haar fluctuation spectrum
#'
#' @description Haar fluctuation spectrum after Lovejoy
#' still needs to be documented
#' see also parcival and walden wavelets book
#' currently requires n=2^i data points
#' @param x the input series (typically a numeric)
#' @param q the fluctuation order
#' @param discarded_scales the number of discarded scales on the lowest frequency side of the dyadic scale
#' @importFrom stats lm
#' @examples
#' x <- rnorm(2048)
#' xi1 <- haar(x, q = 1)
#' xi2 <- haar(x, q = 2)
#' @export
haar <- function(x, q = 2, discarded_scales = 4) {
  n <- length(x)
  nsteps <- floor(log(n) / log(2))
  nmax <- 2^nsteps
  DT <- sapply(seq(nsteps), function(i) {
    m <- 2^i
    nn <- nmax / m
    # organises the series into m columns
    # e.g. for scale 4, we will do 8 columns
    M <- matrix(x, nmax / m, m, byrow = TRUE)
    # take the average over the first m/2 columns
    M1 <- M[, 1:(m / 2), drop = FALSE]
    tM1 <- apply(M1, 1, mean)
    # take the average over last m/2 columns
    M2 <- M[, (m / 2 + 1):m, drop = FALSE]
    tM2 <- apply(M2, 1, mean)
    # take the difference : this effectively
    # computes the convolution by a Haar wavelet of 
    # scale m/2, over successive steps spaced by m
    abs(tM2 - tM1)
  })
  x <- 2^(seq(nsteps)-1)
  # y is the _mean_ convolution of the signal by the Haar 
  # wavelet
  y <- sapply(DT, function(x) mean(x^q))
  H <- data.frame(x = x, y = y, logx = log(x), logy = log(y))
  kept_scales <- seq(max(2, nsteps - discarded_scales))
  out <- lm(H$logy[kept_scales] ~ H$logx[kept_scales])$coefficients[2]
  attr(out, "scale") <- H
  attr(out, "q") <- q
  attr(out, "class") <- "fluctuation_spectrum"
  out
}

#' @rdname haar
#' @param ... Further parameters passed to the plot function
#' @export
plot.fluctuation_specturm <- function(x,...) {
  plot(attr(x, "scale")$x, attr(x, "scale")$y, log = "xy", xlab = "scale", ylab = "Amplitude")
}
