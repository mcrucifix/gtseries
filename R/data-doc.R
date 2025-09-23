#' harmonic_sample
#'
#' A simple harmonic time series with multiple frequency components.
#'
#' @format A list of class \code{discreteSpectrum} with the following components:
#' \describe{
#'   \item{Amp}{Vector of amplitudes of the harmonic components}
#'   \item{Freq}{Vector of frequencies of the harmonic components}
#'   \item{Phases}{Vector of phases of the harmonic components}
#' }
#' The raw time series data are stored as an attribute "data" (a \code{ts} object).
"harmonic_sample"

#' harmonic_sample_noisy
#' 
#' A noisy version of `harmonic_sample`.
#'
#' @format Same as \code{harmonic_sample}, with added Gaussian noise in the "data" attribute.
"harmonic_sample_noisy"

#' pseudo_log
#'
#' A pseudo-log dataset simulating a simple lithological record.
#'
#' @format A list of class \code{discreteSpectrum} with amplitudes, frequencies, and phases.
#' The raw binary time series is stored in the "data" attribute.
"pseudo_log"

#' insol_scaled
#'
#' Scaled test insolation dataset used for time series analysis.
#'
#' @format A \code{ts} object of length 4096.
"insol_scaled"

#' insol_truncated
#'
#' Truncated/thresholded version of \code{insol_scaled}, values in \{-1, 0, 1\}.
#'
#' @format A \code{ts} object of length 4096.
"insol_truncated"
