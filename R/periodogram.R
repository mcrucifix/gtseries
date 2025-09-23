#' Periodogram
#'
#' Simple periodogram
#' @param x : numeric vector or time series object
#' @importFrom stats start
#' @importFrom stats as.ts
#' @author Michel Crucifix for the R code
#

#' @export periodogram
periodogram <- function(x){
  x = stats::as.ts(x)
  dt = deltat(x)
  startx = stats::start(x)[1]
  N <- length(x)
  N2 <- ceiling(N/2)
  freqs <- ((seq(N)-1)/dt/N)[0:N2]
  f <-  fft(x)[0:N2]
  Power <- Mod(f)^2 / (dt*N2)^2
  Phase <- Arg(f)
  OUT <- list(Freq=freqs, Power=Power, Phase=Phase)
  attr(OUT, "class") = "periodogram"
  return(OUT)
}

#' @export periodogram_complex
periodogram_complex <- function(x){
  x = stats::as.ts(x)
  dt = deltat(x)
  startx = stats::start(x)[1]
  N <- length(x)
  N1 <- N-1
  freqs <- ((seq(N)-1)/dt/N)[0:N1]
  f <-  fft(x)[0:N1]

  negatives <- which(freqs > 1/(2*dt))
  freqs[negatives] = freqs[negatives] - 1./dt
  O <- order(freqs)
  freqs <- freqs[O]
  f <- f[O]
   
  Power <- Mod(f)^2 / (dt)^2
  Phase <- Arg(f)
  OUT <- list(Freq=freqs, Power=Power, Phase=Phase)
  attr(OUT, "class") = "periodogram_complex"
  OUT
}


#' Plot a periodogram object
#'
#' @param x periodogram object
#' @param ... further arguments passed to the default plot method, may include:
#'   \describe{
#'     \item{log}{log plot.  axis to be logged, e.g. \code{"x"}, \code{"y"}, or \code{"xy"} (default: \code{"xy"})}
#'     \item{xlab}{x-axis label (default: \code{"Frequency"})}
#'     \item{ylab}{y-axis label (default: \code{"Power density"})}
#'     \item{type}{plot type, e.g. \code{"l"} for lines}
#'     \item{main}{plot title}
#'   }
#' @return The plot
#' @rdname periodogram
#' @method plot periodogram
#' @export
plot.periodogram <- function(x, ...) {
  args <- list(...)
  log   <- if (!is.null(args$log)) args$log else "xy"
  xlab  <- if (!is.null(args$xlab)) args$xlab else "Frequency"
  ylab  <- if (!is.null(args$ylab)) args$ylab else "Power density"
  type  <- if (!is.null(args$type)) args$type else "l"
  plot(x$Freq, x$Power, type=type, log=log, xlab=xlab, ylab=ylab, ...)
}

#' Plot a complex periodogram object
#'
#' @param x periodogram_complex object
#' @param ... further arguments passed to the default plot method, may include:
#'   \describe{
#'     \item{log}{log plot, e.g. \code{"x"}, \code{"y"}, or \code{"xy"} (default: \code{"y"})}
#'     \item{xlab}{x-axis label (default: \code{"Frequency"})}
#'     \item{ylab}{y-axis label (default: \code{"Power density"})}
#'     \item{type}{plot type, e.g. \code{"l"} for lines}
#'     \item{main}{plot title}
#'   }
#' @return The plot
#' @rdname periodogram_complex
#' @method plot periodogram_complex
#' @export
plot.periodogram_complex <- function(x, ...) {
  args <- list(...)
  log   <- if (!is.null(args$log)) args$log else "y"
  xlab  <- if (!is.null(args$xlab)) args$xlab else "Frequency"
  ylab  <- if (!is.null(args$ylab)) args$ylab else "Power density"
  type  <- if (!is.null(args$type)) args$type else "l"
  plot(x$Freq, x$Power, type=type, log=log, xlab=xlab, ylab=ylab, ...)
}

