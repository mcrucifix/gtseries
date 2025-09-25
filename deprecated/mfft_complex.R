#' Modified Fourier transform
#'
#' Implementation of the Frequency Modified Fourier Transform
#' (Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137).
#' Given a quasi-periodic complex signal X + iY, the algorithm
#' estimates the frequencies (f_j), amplitudes (A_j) and phases
#' (psi_j) in its decomposition. The algorithm was generalised
#' by M. Crucifix when the input signal is real.
#' `mfft` is a wrapper that redirects the data
#' to `mfft_real` or `mfft_complex` depending on the nature of
#' the input data.
#' TODO: add Laskar reference and say that the correction=0 version
#'       corresponds to Laskar's routine.
#' @param xdata The data provided either as a time series (advised), or as a vector.
#' may be complex
#' @param min_freq,max_freq If provided, bracket the frequencies to be probed. Note these are
#'        more precisely angular velocities (\eqn{2\pi/\mathrm{period}}), expressed in time-inverse units
#'        with the time resolution encoded in `xdata` if the latter is a time series.
#'        The computed frequencies are in the range given by min_freq and max_freq.
#' @param correction 0: no frequency correction (equivalent to Laskar); 1 : frequency correction using linear approximation ; 2: frequency correction using synthetic data;
#' 3: second order-correction using synthetic data (all documented in the Sidlichovsky and Nesvorny reference). Note: 1 is not implemented for complex time series; 4 is
#' not documented for real time series
#' @param nfreq The number of frequencies returned, must be smaller than the length of xdata.
#' @param force_complex Use the complex number implementation even if the time series is real.
#' @return A `discreteSpectrum` object, based on a data.frame with columns "Freq", "Ampl" and "Phases".
#'         Note that because of a language glitch (to be fixed), "Freq" actually means "Rate".
#' @author Michel Crucifix for the R code, and David Nesvorny for most of the supporting C code doing the
#'         Frequency Modified Fourier transform for Complex Numbers
#' @references Sidlichovsky, M. & Nesvorny, D. (1997). Frequency Modified Fourier Transform. Celestial Mechanics and Dynamical Astronomy, 65, 137.
#' @export mfft
mfft <- function(xdata, nfreq = 15, min_freq = NULL, max_freq = NULL, correction = 1, force_complex = FALSE) {
  if (is.complex(xdata) || force_complex) {
    return(mfft_complex(xdata, nfreq, min_freq, max_freq, correction))
  } else {
    return(mfft_real(xdata, nfreq, min_freq, max_freq, correction))
  }
}

#' Frequency Modified Fourier transform for Complex Numbers
#'
#' Implementation of the Frequency Modified Fourier Transform
#' (Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137).
#' Given a quasi-periodic complex signal X + iY, the algorithm
#' estimates the frequencies (f_j), amplitudes (A_j) and phases
#' (psi_j) in its decomposition.
#' @useDynLib gtseries
#' @importFrom stats start
#' @importFrom stats as.ts
#' @param xdata The data provided either as a time series (advised), or as a vector.
#' may be complex
#' @param min_freq,max_freq If provided, bracket the frequencies to be probed. Note these are
#'        more precisely angular velocities (\eqn{2\pi/\mathrm{period}}), expressed in time-inverse units
#'        with the time resolution encoded in `xdata` if the latter is a time series.
#' @param correction
#'      Modified Fourier Transform                  if   correction = 0;
#'      Frequency Modified Fourier Transform        if   correction = 1;
#'      FMFT with additional non-linear correction  if   correction = 2
#' (while the first algorithm is approximately 3 times faster than the third one,
#' the third algorithm should be in general much more precise).
#' The computed frequencies are in the range given by min_freq and max_freq.
#' @param nfreq The number of frequencies returned, must be smaller than the length of xdata.
#' @return A `discreteSpectrum` object, based on a data.frame with columns "Freq", "Ampl" and "Phases".
#'         Note that because of a language glitch (to be fixed), "Freq" actually means "Rate".
#' @author Michel Crucifix for the R code, and David Nesvorny for most of the supporting C code doing the
#' actual computations
#' @references Sidlichovsky, M. & Nesvorny, D. (1997). Frequency Modified Fourier Transform. Celestial Mechanics and Dynamical Astronomy, 65, 137.
#' @export mfft_complex
mfft_complex <- function(xdata, nfreq = 30, min_freq = NULL, max_freq = NULL, correction = 1) {
  ts_data <- stats::as.ts(xdata)
  time_step <- deltat(ts_data)
  time_start <- stats::start(ts_data)[1]
  imag_data <- Im(ts_data)
  real_data <- Re(ts_data)
  data_length <- length(real_data)

  correction_flag <- c(1, 0, 2, 3)[correction + 1]
  print("correction_flag")
  print(correction_flag)
  if (correction_flag == 0) {
    message("This correction scheme is not implemented for complex time series; will use synthetic data correction instead.")
    correction_flag <- 2
  }

  if (is.null(min_freq)) {
    min_freq <- -pi
    max_freq <- pi
  } else {
    min_freq <- min_freq * time_step
    max_freq <- max_freq * time_step
  }

  spectrum_freq <- as.integer(nfreq)
  spectrum_min_freq <- as.double(min_freq)
  spectrum_max_freq <- as.double(max_freq)
  spectrum_real <- as.double(real_data)
  spectrum_imag <- as.double(imag_data)
  spectrum_signal1 <- as.double(rep(0, 3 * nfreq))
  spectrum_signal2 <- as.double(rep(0, 3 * nfreq))
  spectrum_signal3 <- as.double(rep(0, 3 * nfreq))

  result <- .C(
    "fmft",
    spectrum_freq,
    spectrum_min_freq,
    spectrum_max_freq,
    as.integer(correction_flag),
    as.integer(data_length),
    spectrum_real,
    spectrum_imag,
    spectrum_signal1,
    spectrum_signal2,
    spectrum_signal3,
    DUP = TRUE
  )

  result_matrix <- t(matrix(result[[7 + correction_flag]], 3, nfreq))

  freq_values <- result_matrix[seq(nfreq)] / time_step
  amplitude_values <- result_matrix[nfreq + seq(nfreq)]
  phase_values <- result_matrix[2 * nfreq + seq(nfreq)] - time_start * freq_values

  spectrum <- data.frame(Freq = freq_values, Amp = amplitude_values, Phases = phase_values)
  class(spectrum) <- c("discreteSpectrum", "data.frame")
  attr(spectrum, "nfreq") <- nfreq
  attr(spectrum, "data") <- real_data

  return(spectrum)
}
