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
  Power <- Mod(f)^2 / (dt*N2)^2
  Phase <- Arg(f)
  OUT <- list(Freq=freqs, Power=Power, Phase=Phase)
  attr(OUT, "class") = "periodogram"
  return(OUT)
}

#' @export periodogram_complex
periodogram_complex <- function(xdata){
  xdata = stats::as.ts(xdata)
  dt = deltat(xdata)
  startx = stats::start(xdata)[1]
  N <- length(xdata)
  N1 <- N-1
  freqs <- ((seq(N)-1)/dt/N)[0:N1]
  f <-  fft(xdata)[0:N1]

  negatives <- which(freqs > 1/(2*dt))
  freqs[negatives] = freqs[negatives] - 1./dt
  O <- order(freqs)
  freqs <- freqs[O]
  f <- f[O]
   
  Power <- Mod(f)^2 / (dt)^2
  Phase <- Arg(f)
  OUT <- list(Freq=freqs, Power=Power, Phase=Phase)
  attr(OUT, "class") = "periodogram_complex"
  return(OUT)
}



#' @rdname periodogram
#' @export
plot.periodogram <- function(X,log='xy', xlabel= "Frequency", ylabel="Power density", ...){
  plot(X$Freq, X$Power, type='l', log=log,xlabel=xlabel, ylabel=ylabel, ...)
}

#' @rdname periodogram_complex
#' @export
plot.periodogram_complex <- function(X,log='y', xlabel= "Frequency", ylabel="Power density", ...){
  plot(X$Freq, X$Power, type='l', log=log,xlabel=xlabel, ylabel=ylabel, ...)
}

