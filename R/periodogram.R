#' Periodogram (for real and complex time series)
#'
#' Computes the periodogram for a time series. The default \code{periodogram}
#' is for real-valued time series (returns only non-negative frequencies).
#' The alias \code{periodogram_complex} is provided for complex-valued time
#' series, and returns both positive and negative frequencies.
#'
#' @param x Numeric vector or time series object. For \code{periodogram_complex},
#' may be complex-valued.
#' @return A list containing:
#'   \describe{
#'     \item{Freq}{Frequencies. For \code{periodogram}, non-negative; for
#'       \code{periodogram_complex}, both negative and positive.}
#'     \item{Power}{Power spectral density.}
#'     \item{Phase}{Phase spectrum.}
#'   }
#' @details
#' \code{periodogram} computes the simple periodogram for real-valued signals,
#' returning only the positive frequencies.
#'
#' \code{periodogram_complex} is an alias for use with complex-valued signals.
#' It returns both positive and negative frequencies, as is standard for
#' Fourier transforms of complex time series.
#'
#' @aliases periodogram periodogram_complex
#' @author Michel Crucifix
#' @seealso \code{\link{plot.periodogram}}, \code{\link{plot.periodogram_complex}}
#' @examples
#' # Real-valued signal
#' x <- rnorm(128)
#' pg <- periodogram(x)
#' plot(pg)
#'
#' # Complex-valued signal
#' y <- rnorm(128) + 1i * rnorm(128)
#' pgc <- periodogram_complex(y)
#' plot(pgc)
#'
#' @export periodogram
#' @export periodogram_complex
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

