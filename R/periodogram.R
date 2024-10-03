#' Periodogram
#'
#' Simple periodogram
#' @param xdata 
#' @importFrom stats start
#' @importFrom stats as.ts
#' @author Michel Crucifix for the R code
#

#' @export periodogram
periodogram <- function(xdata){
  xdata = stats::as.ts(xdata)
  dt = deltat(xdata)
  startx = stats::start(xdata)[1]
  N <- length(xdata)
  N2 <- ceiling(N/2)
  freqs <- ((seq(N)-1)/dt/N)[0:N2]
  f <-  fft(xdata)[0:N2]
  Power <- Mod(f)^2
  Phase <- Arg(f)
  out <- list(Freq=freqs, Power=Power, Phase=Phase)
  attr(OUT, "class") = "periodogram"
  return(OUT)
}

#' @rdname periodogram
#' @export
plot.periodogram <- function(X,...){
  plot(X$Freq, X$Mod, type='l', log='xy',...)
}

