# TO DO: THERE IS A BIAS FOR THE LONG LENGHTSCALES
# I GUESS MUST MULTIPLY BY (N-1)/N where N is
# length / scale



#' Haar fluctuation spectrum
#'
#' @description Haar fluctuation spectrum after Lovejoy
#' still needs to be documented
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
    M <- matrix(x, nmax / m, m, byrow = TRUE)
    M1 <- M[, 1:(m / 2), drop = FALSE]
    M2 <- M[, (m / 2 + 1):m, drop = FALSE]
    tM1 <- apply(M1, 1, mean)
    tM2 <- apply(M2, 1, mean)
    abs(tM2 - tM1) * (2 * nn) / (2 * nn - 1)
  })
  x <- 2^seq(nsteps)
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
