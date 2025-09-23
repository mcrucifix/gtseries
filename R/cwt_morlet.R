#' Continuous Morlet Wavelet Transform and Related Methods
#'
#' This file provides functions for performing the Continuous Morlet Wavelet Transform (CWT)
#' on input time series, visualizing the results, extracting ridges, and reconstructing time series.
#'
#' References:
#' Mallat, S. 1998: A wavelet Tour of Signal Processing. Academic Press, 577 pp.
#' Tchawitchia P., Wavelets, functions and operators, ch. 3 in Wavelets : theory and applications,
#' Erlenbacher et al. eds, Oxford University Press 1996 [for ridge extraction]
#'
#' @importFrom stats fft deltat time ts qnorm
#' @importFrom graphics axis image lines par text
#' @importFrom fields image.plot
#' @importFrom signal unwrap
#' @keywords internal

col_wavelet <- colorRampPalette(c("darkblue", "lightblue", "grey", "white", "red"))(100)

#' Ridge Extraction from CWT
#'
#' Extracts the ridge (maximal amplitude curve) from a continuous wavelet transform.
#'
#' @param A Numeric or time series vector.
#' @param first_guess Initial scale guess.
#' @param k0 Morlet parameter.
#' @param tol Tolerance for convergence.
#' @param itmax Maximum iterations.
#' @param B Time index vector.
#' @param window Numeric vector for scale window.
#' @return Ridge amplitude vector with attributes.
search_ridge <- function(A, first_guess, k0 = 5.6, tol = 1e-6, itmax = 20, B = seq(along = A), window = c(0, Inf)) {
  n <- length(A)
  amp <- A
  amp[] <- NA
  ome <- amp
  scale_vec <- amp
  dt <- deltat(A)
  scale_val <- first_guess
  for (b in B)
  {
    delta <- Inf
    iter <- 0
    while ((delta > tol) & (iter <= itmax)) {
      w <- cwt_morlet(A, k0 = k0, calcmask = FALSE, scale = scale_val, deriv = TRUE)
      p <- attr(w, "deriv")[b]
      expected_scale <- k0 / Re(p)
      delta <- scale_val - expected_scale
      scale_val <- expected_scale
      iter <- iter + 1
      if (is.na(scale_val)) {
        iter <- itmax
      }
    }
    if (iter < itmax & scale_val > window[1] & scale_val < window[2]) {
      amp[b] <- w[b]
      ome[b] <- p
      scale_vec[b] <- scale_val
    } else {
      amp[b] <- NA
      ome[b] <- NA
      scale_vec[b] <- NA
      if (scale_val > window[2]) scale_val <- window[2]
      if (scale_val < window[1]) scale_val <- window[1]
    }
  }
  attr(amp, "period") <- 2 * pi * dt / ome
  attr(amp, "scale") <- scale_vec
  attr(amp, "k0") <- k0
  attr(amp, "class") <- "ridge"
  amp
}

#' Enhance Ridge
#'
#' Enhances a previously extracted ridge.
#'
#' @param A Input time series.
#' @param R Ridge object.
#' @param confidence Confidence interval for enhancement.
#' @param inter Step size for enhancement.
#' @return Enhanced ridge amplitude vector.
enhance_ridge <- function(A, R, confidence = 0.95, inter = 0.08) {
  scale_vec <- attr(R, "scale")
  amp <- R
  period_mat <- matrix(NA, nrow = length(A), ncol = 3)
  k0 <- attr(R, "k0")
  dt <- deltat(A)

  width <- log(qnorm(0.5 + confidence / 2) * sqrt(2) / k0 + 1)

  for (b in seq(along = R)) {
    print(b)
    if (!is.na(scale_vec[b])) {
      scale_seq <- exp(log(scale_vec[b]) + seq(-width * 1.03214, +width, inter))
      ## 1.03214 is an empirical factor to take into account the ridge asymmetry
      w <- cwt_morlet(A, k0 = k0, calcmask = FALSE, scale = scale_seq, deriv = TRUE)
      amp[b] <- sum(w[b, ]) * k0 / sqrt(2 * pi) * inter
      period_mat[b, ] <- 2 * pi * dt / k0 * c(scale_vec[b], min(scale_seq), max(scale_seq))
    }
  }
  attr(amp, "period") <- period_mat
  amp
}

#' Reconstruct Time Series from Morlet Wavelet Transform
#'
#' Reconstructs the original time series from its Morlet wavelet transform, optionally
#' restricting to certain scales or periods.
#'
#' @param W A wavelet object.
#' @param scales Numeric vector. Restrict reconstruction to these scales.
#' @param periods Numeric vector. Restrict to these periods.
#' @return Time series object.
#' @export
#' @importFrom stats ts
reconstruct_morlet <- function(W, scales = c(-Inf, Inf), periods = NULL) {
  if (!(attr(W, "class") == "wavelet")) stop("object is not a wavelet transform")
  if (!(attr(W, "wavelet") == "morlet")) stop("object is not a MORLET  wavelet transform")
  a <- attr(W, "scale")
  j <- which(a > scales[1] & a < scales[2])
  # if period is set, then overrides scales
  if (!is.null(periods)) {
    a <- attr(W, "period")
    j <- which(a > periods[1] & a < periods[2])
  }

  inter <- attr(W, "parameters")$inter
  k0 <- attr(W, "parameters")$k0

  output <- rowSums(W[, j]) * log(2) / inter / (sqrt(2 * pi)) * k0
  output <- ts(output, start = attr(W, "time")[1], deltat = diff(attr(W, "time"))[1])
}

#' Continuous Cross Morlet Wavelet Transform
#'
#' Computes the cross wavelet transform between two time series.
#'
#' @param A First time series.
#' @param B Second time series.
#' @param ... Additional parameters for Morlet transform.
#' @return Cross wavelet transform matrix.
cross_morlet <- function(A, B, ...) {
  CA <- cwt_morlet(A, ...)
  CB <- cwt_morlet(B, ...)
  (CA * Conj(CB)) / (Mod(CA) * Mod(CB))
}

#' Continuous Morlet Wavelet Transform
#'
#' Computes the Continuous Morlet Wavelet Transform (CWT) of a time series.
#'
#' Returns a wavelet object containing the complex CWT coefficients, with attributes for time, scale, period, and more.
#'
#' @param A Numeric or time series vector.
#' @param inter Integer. Number of intervals per octave.
#' @param k0 Numeric. Morlet wavelet parameter.
#' @param amin Minimum scale.
#' @param amax Maximum scale.
#' @param calcmask Logical. Whether to mask edge effects.
#' @param scale Numeric. Specific scale(s) to use.
#' @param deriv Logical. Whether to calculate derivatives.
#' @return A wavelet object (matrix with attributes).
#' @importFrom stats fft deltat time
#' @importFrom graphics image axis par
#' @export
cwt_morlet <- function(A, inter = 20, k0 = 5.6, amin = 1, amax = Inf, calcmask = TRUE, scale = NA, deriv = FALSE) {
  y <- A
  time_vec <- stats::time(A)
  dt <- deltat(A)
  n_y <- length(y)

  if (is.na(scale)) {
    local({
      n_y2 <- round(n_y / 2)
      exp1 <- log2(amin) + 2
      exp2 <- min(round(log2(n_y2)) + 1, amax)

      scale <- vector()

      j <- 0
      for (m in seq(exp1, exp2 - 1)) {
        jj <- inter - 1
        for (n in seq(0, jj)) {
          a <- 2^(m + n / inter)
          j <- j + 1
          scale[j] <<- a
        }
      }
    })
  }

  omega0 <- k0 / scale
  period <- 1 / omega0 * 2 * pi * dt

  x <- y
  n <- length(x)

  k <- (0:(n - 1)) * 2 * pi / n

  f <- stats::fft(x)

  J <- length(scale)

  wave <- matrix(as.complex(0), nrow = n, ncol = J)

  if (calcmask) mask <- matrix(TRUE, nrow = n, ncol = J)
  if (deriv) dwave <- wave

  # n_k <- length(k)
  for (a1 in seq(along = scale)) {
    expnt <- -(scale[a1] * k - k0)^2 / 2
    norm <- 2
    daughter <- norm * exp(expnt)
    wave[, a1] <- stats::fft(f * daughter, inverse = TRUE) / n
    if (deriv) dwave[, a1] <- stats::fft(f * daughter * (-(1i * k)), inverse = TRUE) / n
    if (calcmask) {
      mask[1:min(n, ceiling(sqrt(2) * scale[a1])), a1] <- NA
      mask[(n - min(n, ceiling(sqrt(2) * scale[a1]))):n, a1] <- NA
    }
  }
  if (deriv) {
    dwave <- dwave / (-1i * wave)
  }
  time_vec <- as.array(time_vec)
  period <- as.array(period)
  attr(wave, "time") <- time_vec
  if (calcmask) attr(wave, "mask") <- mask
  attr(wave, "period") <- period
  attr(wave, "scale") <- scale
  attr(wave, "class") <- "wavelet"
  attr(wave, "wavelet") <- "morlet"
  attr(wave, "parameters") <- list(k0 = k0, inter = inter, deltat = dt)
  if (deriv) attr(wave, "deriv") <- dwave
  wave
}

#' Plot a Wavelet Object
#'
#' Visualizes the modulus (amplitude) of a wavelet transform in time-period space.
#'
#' @param x A wavelet object as returned by \code{cwt_morlet}.
#' @param resx Resolution in time axis.
#' @param resy Resolution in period axis.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param scaling_correction Numeric. Power scaling correction factor.
#' @param col Color palette.
#' @param legend Logical. Whether to include a legend.
#' @param axes Logical. Whether to plot axes.
#' @param Mode Function to compute modulus.
#' @param plotMask Logical. Whether to plot cone of influence mask.
#' @param ... Additional arguments passed to image plotting functions.
#' @importFrom fields image.plot
#' @importFrom graphics image axis par
#' @export
plot.wavelet <- function(x, resx = 400, resy = 300,
                         xlab = "Time", ylab = "Period", scaling_correction = 0,
                         col = col_wavelet, legend = FALSE,
                         axes = TRUE, Mode = Mod, plotMask = TRUE, ...) {
  require(fields)
  time_vec <- attr(x, "time")
  period <- attr(x, "period")
  mask <- attr(x, "mask")
  if (!plotMask) mask[, ] <- 1 # do not hide influence cone
  scale_vec <- attr(x, "scale")
  thin_factor_x <- max(ceiling(length(time_vec) / resx), 1)
  thin_factor_y <- max(ceiling(length(scale_vec) / resy), 1)
  sub_x <- seq(thin_factor_x, length(time_vec), thin_factor_x)
  sub_a <- seq(thin_factor_y, length(scale_vec), thin_factor_y)
  wave_scaled <- Mode(x[sub_x, sub_a]) * mask[sub_x, sub_a]
  for (i in 1:length(sub_a)) wave_scaled[, i] <- wave_scaled[, i] / (scale_vec[i]^scaling_correction)
  if (legend) par(oma = c(2, 2, 2, 5))
  image(time_vec[sub_x], period[sub_a], wave_scaled, log = "y", ylab = ylab, xlab = xlab, axes = FALSE, col = col, ...)
  if (axes) {
    axis(1)
    axis.log10(2, "")
  }
  if (legend) {
    par(oma = c(2, 2, 2, 2))
    image.plot(time_vec[sub_x], period[sub_a], wave_scaled,
               legend.only = TRUE, log = "y",
               ylab = ylab, xlab = xlab,
               axes = FALSE, col = col, ...)
    par(oma = c(2, 2, 2, 5))
  }
}

#' Power Spectrum of a Wavelet Transform
#'
#' Computes the mean power spectrum from a wavelet transform.
#'
#' @param wave A wavelet object.
#' @return Data frame with frequency and power columns.
#' @export
powerspectrum.wavelet <- function(wave, ...) {
  scale_vec <- attr(wave, "scale")
  period <- attr(wave, "period")
  time_vec <- attr(wave, "time")
  freq <- 1. / period
  power <- vector("double", length(freq))
  for (j in seq_along(freq)) power[j] <- mean(Mod(wave[, j])^2, na.rm = TRUE)
  data.frame("frequency" = freq, "power" = power)
}

#' Analytic Ridge Extraction
#'
#' Performs analytic ridge extraction from a wavelet transform.
#'
#' @param wave Wavelet object.
#' @param ridge_only Logical. Return only ridge if TRUE.
#' @param plim Numeric vector. Period limits.
#' @return Ridge or modified wavelet object.
#' @importFrom signal unwrap
analytic.ridge <- function(wave, ridge_only = FALSE, plim = c(-Inf, +Inf)) {
  require("signal")
  period <- attr(wave, "period")
  n_period <- length(period)
  period_out <- which(period < plim[1] | period > plim[2])
  scale_vec <- attr(wave, "scale")
  time_vec <- attr(wave, "time")
  dt <- diff(attr(wave, "time"))[1]
  WR <- wave
  WR[, period_out] <- NA
  # scale for scalogram
  for (i in seq(1, n_period)) WR[, i] <- WR[, i] / sqrt(scale_vec[i])

  WP <- Arg(WR)
  for (i in seq(1, n_period)) WP[, i] <- unwrap(WP[, i])
  for (i in seq(1, n_period)) WP[, i] <- c(diff(WP[, i]), NA) / (2 * pi) / dt * period[i]
  WR[which(abs(WP - 1) > 0.025)] <- NA
  if (ridge_only) {
    R <- vector()
    for (i in seq(1, length(time_vec))) {
      M <- which.max(Mod(WR[i, ]))
      if (any(M)) R[i] <- WR[i, M]
    }
    R
  } else {
    WR
  }
}
