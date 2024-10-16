#' harmonic_sample
#'
#' A simple time series for testing purposes
#' cos(t*0.13423167+0.00) + 1.3 * cos(t*0.119432+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) 
#' for t from 0 to 1023. It is coded as a `discreteSpectrum` objcet and self-describes its spectral decomposition. 
#'
#' @docType data
#' @keywords datasets
#' @name harmonic_sample
#' @usage data(harmonic_sample)
NULL

#' harmonic_sample_noisy
#'
#' Same has harmonic sample but with an added random time series
#' cos(t*0.13423167+0.00) + 1.3 * cos(t*0.119432+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) + rnorm(1024)*0.12
#' for t from 0 to 1023. It is coded as a `discreteSpectrum` objcet and self-describes its spectral decomposition. 
#'
#' @docType data
#' @keywords datasets
#' @name harmonic_sample_noisy
#' @usage data(harmonic_sample_noisy)
NULL

#' pseudo_log
#'
#' A truncated version (squared-shape) version of a harmonic signal for mimicking a sedimentary log
#' sign( cos(t*0.13423167+0.00) + 1.3 * cos(t*0.119432+2.314) + 0.134994 + 0.4*cos(t*0.653167) + 0.11 * cos(t*0.78913498) + 0.8))
#' for t from 0 to 8191. It is coded as a `discreteSpectrum` objcet and self-describes its spectral decomposition. 
#' @docType data
#' @keywords datasets
#' @name pseudo_log
#' @usage data(pseudo_log)
NULL
