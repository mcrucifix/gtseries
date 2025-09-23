#' Develop a spectrum into a time series (generic)
#'
#' This is a generic function to reconstruct a time series from a spectrum decomposition.
#'
#' @param M input object (typically a spectrum or decomposition)
#' @param start If supplied, the start time of the decomposition (overrides `times`)
#' @param end If supplied, the end time of the decomposition (used with `start` and `deltat`)
#' @param deltat If supplied, the time interval (used with `start` and `end`)
#' @param times If supplied, explicit times for the decomposition
#' @param ... Further arguments passed to methods
#' @note Place holder for type-specific develop functions.
#' @return Nothing; method dispatches to appropriate method.
#' @export
develop <- function(M, start = NULL, end = NULL, deltat = NULL, times = NULL, ...) {
  UseMethod("develop")
}

#' Trigonometric exponential function
#'
#' Returns the complex exponential of x.
#'
#' @param x Numeric input value.
#' @return Complex exponential of x.
#' @export
cis <- function(x) exp(1i * x)

#' Discrete spectrum reconstruction
#'
#' Reconstructs a time series (or its components) from a discrete spectrum decomposition.
#'
#' @param M A `discreteSpectrum` object
#' @param start If supplied, the start time of the decomposition (overrides `times`)
#' @param end If supplied, the end time of the decomposition (used with `start` and `deltat`)
#' @param deltat If supplied, the time interval (used with `start` and `end`)
#' @param times If supplied, explicit times for the decomposition
#' @param dfunction Trigonometrical function, e.g. 'cos', 'sin', or 'cis'
#' @param maxfreq Maximum frequency to include in the reconstruction
#' @param sum TRUE if the user wants to sum components in the reconstruction
#' @param trendshift TRUE if the user wants to account for trend and shift encoded in the `discreteSpectrum` object
#' @note If none of `times`, `start`, and `deltat` are supplied, will reconstruct based on the attribute `data` (which must be present). If no `data` attribute is available, an error is returned.
#' @return If `sum = FALSE`, a list of reconstructed components. If `sum = TRUE`, the full reconstructed time series.
#' @export
develop.discreteSpectrum <- function(M, start = NULL, end = NULL, deltat = NULL, times = NULL,
                                    dfunction = cos, maxfreq = NULL, sum = TRUE, trendshift = TRUE) {
  if (!("discreteSpectrum" %in% class(M))) stop("object is not a discreteSpectrum decomposition")

  times_is_ts <- FALSE
  if (is.ts(times)) {
    start <- start(times)
    deltat <- deltat(times)
    end <- start + length(times) * deltat
    times_is_ts <- TRUE
  }
  if (!is.null(start)) {
    if (is.null(deltat) || is.null(end)) stop("if you supply start, you must also supply deltat and end")
    n <- (end - start) %/% deltat
    times <- ts(start + seq(0, n) * deltat, start = start, deltat = deltat)
    times_is_ts <- TRUE
  }

  if (is.null(times)) {
    if (is.null(attr(M, "data"))) stop("if you do not supply any time argument (times, or (start, end, deltat)), then object must have a valid data attribute")
    x_data <- attr(M, "data")
    start <- stats::start(x_data)[1]
    deltat <- stats::deltat(x_data)
    times <- (seq(length(x_data)) - 1) * deltat + start
    times_is_ts <- TRUE
  }

  n_freq <- attr(M, "nfreq")
  if (is.null(n_freq)) n_freq <- length(M$Amp)
  if (!is.null(maxfreq)) n_freq <- min(n_freq, maxfreq)

  if (times_is_ts) {
    reconstructed <- lapply(seq(n_freq), function(i)
      ts(M$Amp[i] * dfunction(M$Freq[i] * times + M$Phase[i]), start = start, deltat = deltat)
    )
  } else {
    reconstructed <- sapply(seq(n_freq), function(i)
      M$Amp[i] * dfunction(M$Freq[i] * times + M$Phase[i])
    )
  }

  if (sum) {
    shift <- attr(M, "shift"); if (is.null(shift)) shift <- 0
    trend <- attr(M, "trend"); if (is.null(trend)) trend <- 0
    if (!trendshift) { shift <- 0; trend <- 0 }
    if (times_is_ts) {
      reconstructed <- Reduce('+', reconstructed) + trend * times + shift
    } else {
      reconstructed <- apply(reconstructed, 1, sum) + trend * times + shift
    }
  } else if (times_is_ts) {
    reconstructed <- lapply(reconstructed, function(x) ts(x, start = start, deltat = deltat))
  }
  return(reconstructed)
}

#' Convert discreteSpectrum object to data.frame
#'
#' Converts a discreteSpectrum object into a data.frame.
#' @param x discreteSpectrum object
#' @return Data.frame with columns Freq, Amp, Phases
#' @export
as.data.frame.discreteSpectrum <- function(x) {
  data.frame(Freq = x$Freq, Amp = x$Amp, Phases = x$Phases)
}

#' Plot discrete spectrum
#'
#' Plots the amplitude spectrum of a discreteSpectrum object.
#' @param M discreteSpectrum object
#' @param periods If TRUE, will add a lower axis with period labels
#' @param labels Labels to be set above the frequency peaks
#' @param ... Further graphical arguments
#' @export
plot.discreteSpectrum <- function(M, periods = FALSE, labels = NULL, ...) {
  plot(abs(M$Freq), abs(M$Amp), 'h', ylab = "Amplitudes", xlab = "", ...)
  if (periods) {
    frequencies <- pretty(range(M$Freq / (2 * pi)))
    plabels <- as.character(1 / frequencies)
    if (0 %in% frequencies) plabels[which(frequencies == 0)] <- "∞"
    axis(1, line = 3, at = 2 * pi * frequencies, labels = plabels)
    mtext("Rate", 1, 2)
    mtext("Period", 1, 4)
  } else {
    mtext("Rate", 1, 3)
  }
  if (!is.null(labels)) {
    y_shift <- 0.05 * diff(range(M$Amp))
    text(M$Freq, M$Amp, labels, srt = 90, adj = -0.4)
  }
}

#' Add lines to discrete spectrum plot
#'
#' Adds lines and points to the spectrum plot for a discreteSpectrum object.
#' @param M discreteSpectrum object
#' @param ... Further graphical arguments
#' @export
lines.discreteSpectrum <- function(M, ...) {
  lines(abs(M$Freq), abs(M$Amp), "h", ...)
  points(abs(M$Freq), abs(M$Amp), "p", ...)
}

#' Print discrete spectrum
#'
#' Prints a summary of the discreteSpectrum object.
#' @param M discreteSpectrum object
#' @param ... Further arguments
#' @export
print.discreteSpectrum <- function(M, ...) {
  n <- nrow(as.data.frame(M))
  to_print <- seq(min(10, n))
  print.data.frame(cbind(as.data.frame(M)[to_print, ], Period = 2 * pi / M$Freq[to_print]))
  if (n > 10) print(sprintf("... + %d other rows \n", n - 10))
}
