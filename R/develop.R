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
#' @importFrom stats is.ts
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
#' Reconstruct a time series from a discreteSpectrum object.
#'
#' @param M discreteSpectrum object
#' @param start see \code{develop}
#' @param end see \code{develop}
#' @param deltat see \code{develop}
#' @param times see \code{develop}
#' @param ... Further arguments passed through. May include:
#'   - `dfunction`: trigonometrical function (`cos`, `sin`, or `cis`)
#'   - `max_freq`: maximum number of frequencies to include
#'   - `sum`: logical; sum components in reconstruction?
#'   - `trendshift`: logical; account for trend and shift encoded in the object?
#' @note If none of times, start and deltat are supplied, will reconstruct based on the attribute `xdata`.
#' @return List of reconstructed components if sum=FALSE, full reconstructed time series otherwise.
#' @method develop discreteSpectrum
#' @export
develop.discreteSpectrum <- function(M, start = NULL, end = NULL, deltat = NULL, times = NULL, ...) {
  if (!("discreteSpectrum" %in% class(M))) stop("object is not a discreteSpectrum decomposition")

  args <- list(...)
  dfunction   <- if (!is.null(args$dfunction)) args$dfunction else cos
  max_freq  <- if (!is.null(args$max_freq)) args$max_freq else NULL
  sum  <- if (!is.null(args$sum)) args$ylab else TRUE
  trendshift  <- if (!is.null(args$trendshift)) args$trendshift else TRUE

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
  max_freq <- attr(M, "max_freq")
  if (!is.null(max_freq)) n_freq <- min(n_freq, max_freq)

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
#' @param row.names kept for compatibility with as.data.frame S3 method
#' @param optional kept for compatibility with as.data.frame S3 method
#' @param ... passed to the data.frame function 
#' @return Data.frame with columns Freq, Amp, Phases
#' @export
as.data.frame.discreteSpectrum <- function(x, row.names, optional, ...) {
  data.frame(Freq = x$Freq, Amp = x$Amp, Phases = x$Phases, ... )
}

#' Plot discrete spectrum
#'
#' Plots the amplitude spectrum of a discreteSpectrum object.
#' @param x discreteSpectrum object
#' @param ... Further arguments passed through. May include:
#'   - `periods`: Add a period axis bar below the frequency axis. Defaults to FALSE
#'   - `labels`: if set, is passed to the standard plot function
#
#' @export
plot.discreteSpectrum <- function(x, ...) {
  args <- list(...)
  periods <- args$periods
  args$periods <- NULL  # Remove 'periods' to avoid sending to plot()

  # ---- Validate 'periods' ----
  if (is.null(periods)) {
    periods <- FALSE
  } else if (!is.logical(periods) || length(periods) != 1) {
    warning("Argument 'periods' must be TRUE or FALSE. Using FALSE.")
    periods <- FALSE
  }

  # ---- Call plot with cleaned args ----
  do.call(
    plot,
    c(
      list(abs(x$Freq), abs(x$Amp), type = 'h', ylab = "Amplitudes", xlab = ""),
      args
    )
  )

  # ---- Custom axis ----
  if (periods) {
    # Detect zero frequency (infinite period case)
    has_zero <- any(abs(x$Freq) < .Machine$double.eps)

    # Compute *finite* periods only
    finite_periods <- 1 / (x$Freq[x$Freq != 0] / (2 * pi))
    periods_pretty <- pretty(range(finite_periods), n = 5)

    # Insert Inf at the beginning if zero frequency exists
    if (has_zero) {
      periods_labels <- c(Inf, periods_pretty)
      freq_at <- c(0, 2 * pi / periods_pretty)
    } else {
      periods_labels <- periods_pretty
      freq_at <- 2 * pi / periods_pretty
    }

    # Convert labels: replace Inf with "∞"
    label_text <- ifelse(is.infinite(periods_labels), "\u221E", as.character(periods_labels))

    axis(1, line = 3, at = freq_at, labels = label_text)
    mtext("Rate", 1, 2)
    mtext("Period", 1, 4)

  } else {
    mtext("Rate", 1, 3)
  }

  # ---- Optional text labels ----
  if (!is.null(args$labels)) {
    # y_shift <- 0.05 * diff(range(x$Amp))
    text(x$Freq, x$Amp, args$labels, srt = 90, adj = -0.4)
  }
}


#' Add lines to discrete spectrum plot
#'
#' Adds lines and points to the spectrum plot for a discreteSpectrum object.
#' @param x discreteSpectrum object
#' @param ... Further graphical arguments
#' @export
lines.discreteSpectrum <- function(x, ...) {
  lines(abs(x$Freq), abs(x$Amp), "h", ...)
  points(abs(x$Freq), abs(x$Amp), "p", ...)
}

#' Print discrete spectrum
#'
#' Prints a summary of the discreteSpectrum object.
#' @param x discreteSpectrum object
#' @param ... Further arguments
#' @export
print.discreteSpectrum <- function(x, ...) {
  n <- nrow(as.data.frame(x))
  to_print <- seq(min(10, n))
  print.data.frame(cbind(as.data.frame(x)[to_print, ], Period = 2 * pi / x$Freq[to_print]))
  if (n > 10) print(sprintf("... + %d other rows \n", n - 10))
}
