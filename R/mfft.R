# Helpers ----------------------------------------------------------------

pisquare <- pi^2
pihalf <- pi / 2

Q <- function(w_T) {
        sin(w_T) / (w_T) * (pisquare) / (pisquare - (w_T * w_T))
}
Q_prime <- function(y) {
        ifelse(y == 0, 0, pisquare / (pisquare - y * y) / y *
               (cos(y) + (sin(y) / y) * (3 * y * y - pisquare) / (pisquare - y * y)))
}
Q_second_0 <- 2 / pisquare - 1. / 3.



cis <- function(x) exp(1i * x)

#' Internal workhorse of the modified Fourier transform for complex and real series
#'
#' Implementation of the Frequency Modified Fourier Transform
#' (Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137).
#' Iteratively extracts dominant frequencies from a real-valued time series
#' using a modified Fourier transform with frequency correction. The function
#' applies a Hanning window to the data, and amplitudes are renormalised through
#' a Gram–Schmidt procedure based on cross-products of Fourier components.
#' The cross-products can be computed either via analytical approximation
#' (faster, following Sidlichovsky and Nesvorny) or explicitly by numerical
#' integration.
#'
#' @param x_data Numeric vector or time series, the real-valued signal to analyse.
#' @param n_freq Integer. Number of frequencies to extract.
#' @param fast Logical (default \code{TRUE}). If \code{TRUE}, use analytical
#'   approximations for Fourier component cross-products; if \code{FALSE},
#'   compute them numerically.
#' @param nu Optional numeric vector of frequencies. If provided, these are used
#'   instead of searching for maxima in the Fourier transform.
#' @param min_freq,max_freq Optional numeric scalars. Real, positive bounds
#'   for the angular frequencies considered in the search. Units depend on
#'   the sampling of \code{x_data}.
#' @return A \code{data.frame} with columns:
#'   \itemize{
#'     \item \code{nu}: estimated angular frequencies
#'     \item \code{amp}: amplitudes of the extracted components
#'     \item \code{phase}: corresponding phases
#'   }
#' The number of rows equals \code{n_freq}.
#'
#' @details This function assumes by design that the input series is real and
#' applies a Hanning window to the data. The iterative Gram–Schmidt
#' renormalisation ensures orthogonality of extracted Fourier components.
#'
#' @references
#' \insertRef{sidlichovsky97aa}{gtseries}
#'
#' @author Michel Crucifix
#' @aliases mfft_real_analyse mfft_complex_analyse
#' @rdname mfft_internal
#' @keywords internal
#
# Workhorse for real series ------------------------------------------------

mfft_real_analyse <- function(x_data, n_freq, fast = TRUE, nu = NULL,
                         min_freq = NULL, max_freq = NULL) {
        if (!is.null(nu)) {
                nu_temp <- unlist(lapply(nu, function(a) if (a == 0) a else c(a, -a)))
                phase <- unlist(lapply(nu, function(a) if (a == 0) 0 else c(0, pihalf)))
                nu <- nu_temp
                if (length(nu) < 2 * n_freq) {
                        nu[2 * n_freq] <- NA
                        phase[2 * n_freq] <- NA
                }
        }   else {
                nu <- rep(NA, 2 * n_freq)
                phase <- rep(NA, 2 * n_freq)
        }

        N <- length(x_data)
        N2 <- N / 2
        hann <- function(N) (1 - cos(2 * pi * seq(0, (N - 1)) / (N)))
        h_N <- hann(N)
        t <- seq(N) - 1
        power <- function(x) (Mod(fft(x))^2)[1:floor((N - 1) / 2)]

        h_prod <- function(x, y) sum(x * h_N * y) / N

        S <- rep(0, 2 * n_freq)
        Prod <- rep(0, 2 * n_freq)
        amp <- rep(NA, 2 * n_freq)
        A <- matrix(0, 2 * n_freq, 2 * n_freq)
        Q_matrix <- matrix(0, 2 * n_freq, 2 * n_freq)
        f <- list()
        x <- list()
        freqs <- 2. * pi * seq(0, (N - 1)) / N
        x[[1]] <- x_data

        for (m in seq(2 * n_freq)) {
                hx <- h_N * x[[m]]
                if (!(m == 2 * n_freq && is.na(nu[m]))) {
                        if (is.na(nu[m])) {
          f_base <- freqs[which.max(power(hx))]
          brackets <- c(f_base - pi / N, f_base + pi / N)
          brackets[1] <- max(min_freq, brackets[1])
          brackets[2] <- min(max_freq, brackets[2])

          tomax <- function(t) {
            function(f) {
              ft <- f * t
              a <- hx %*% cbind(cos(ft), sin(ft))
              a[1] * a[1] + a[2] * a[2]
            }
          }
          f_max <- cmna::goldsectmax(tomax(t),
            brackets[1], brackets[2],
            tol = 1.e-10, m = 9999
          )

         if (f_max > freqs[2] / 2) {
          phase[m] <- 0.
          nu[m] <- f_max
          phase[m + 1] <- pihalf
          nu[m + 1] <- -f_max
        } else {
          phase[m] <- 0.
          nu[m] <- 0.
        }
      }

      f[[m]] <- cos(nu[m] * t + phase[m])

      if (fast) {
        for (i in seq(m)) {
          num <- (nu[m] - nu[i]) * N2
          nup <- (nu[m] + nu[i]) * N2
          phim <- (phase[m] - phase[i])
          phip <- (phase[m] + phase[i])

          Qm <- ifelse(num == 0, 1, Q(num))
          Qp <- ifelse(nup == 0, 1, Q(nup))
          cosm <- cos(phim) * cos(num) * Qm
          sinm <- sin(phim) * sin(num) * Qm
          cosp <- cos(phip) * cos(nup) * Qp
          sinp <- sin(phip) * sin(nup) * Qp

          Q_matrix[m, i] <- 0.5 * (cosm + cosp - sinm - sinp)
        }
      } else {
        for (i in seq(m)) {
          Q_matrix[m, i] <- h_prod(f[[m]], f[[i]])
        }
      }

      A[m, ] <- 0
      A[m, m] <- 1.
      if (m > 1) {
        f_m_bi <- rep(0, (m - 1))
        for (j in seq(m - 1)) for (i in seq(j)) f_m_bi[j] <- f_m_bi[j] - A[j, i] * Q_matrix[m, i]
        # for (j in seq(m - 1)) for (i in seq(j, (m - 1))) A[m, j] <- A[m, j] - f_m_bi[i] * A[i, j]
      }

      # norm <- 0
      # if (m > 1) {
      #   for (i in seq(m)) {
      #     norm <- norm + (A[m, i] * A[m, i] * Q_matrix[i, i])
      #     if (i > 1) {
      #       for (j in seq(i - 1)) {
      #         norm <- norm + 2 * A[m, i] * A[m, j] * Q_matrix[i, j]
      #       }
      #     }
      #   }
      # } else {
      #   norm <- A[m, m] * A[m, m] * Q_matrix[m, m]
      # }
      #

      norm <- Q_matrix[m, m]
      if (m > 1)  for (i in seq(m-1))  norm <- norm - (f_m_bi[i])^2

      # print(sprintf("NORM = %9.4g ", norm - norm2))

      A[m,m] <- A[m,m] / sqrt(norm)

      if (m > 1) {
        for (j in seq(m - 1)) for (i in seq(j, (m - 1))) A[m, j] <- A[m, j] + A[m,m] *  f_m_bi[i] * A[i, j]
      }


      Prod[m] <- h_prod(x[[1]], f[[m]])

      S[m] <- 0.
      for (j in seq(m)) S[m] <- S[m] + Prod[j] * A[m, j]

      x[[m + 1]] <- x[[m]]
      for (j in seq(m)) x[[m + 1]] <- x[[m + 1]] - S[m] * A[m, j] * f[[j]]
    }
  }

  m_max <- 2 * n_freq
  if (is.na(nu[m_max])) m_max <- m_max - 1

  amp[1:m_max] <- 0
  for (m in seq(m_max)) for (j in seq(m)) amp[j] <- amp[j] + S[m] * A[m, j]

  for (m in seq(m_max)) {
    if ((m > 1) && (nu[m - 1] == -nu[m])) {
      phase[m] <- Arg(-1i * amp[m] + amp[m - 1])
      amp[m] <- sqrt(amp[m - 1]^2 + amp[m]^2)
      amp[m - 1] <- NA
      phase[m - 1] <- NA
    }
  }

  valid_frequencies <- which(!is.na(amp))
  nu <- -nu[valid_frequencies]
  amp <- amp[valid_frequencies]
  phase <- phase[valid_frequencies]
  if (length(valid_frequencies) != n_freq) {
          message(sprintf("something goes wrong : %i valid frequencies, and n_freq = %i",
                          valid_frequencies, n_freq))
  }

  OUT <- data.frame(nu = nu, amp = amp, phase = phase)
  return(OUT)
}




# Workhorse for complex series ---------------------------------------------

#' @rdname mfft_internal
#' @keywords internal
mfft_complex_analyse <- function(x_data, n_freq, fast = TRUE, nu = NULL,
                                 min_freq = NULL, max_freq = NULL, use_C_code = TRUE) {
  if (is.null(nu)) {
    nu <- rep(NA, n_freq)
    phase <- rep(NA, n_freq)
  }

  N <- length(x_data)
  N2 <- N / 2
  hann <- function(N) (1 - cos(2 * pi * seq(0, (N - 1)) / (N)))
  h_N <- hann(N)
  t <- seq(N) - 1
  power <- function(x) (Mod(fft(x))^2)

  h_prod <- function(x, y) (x * h_N) %*% Conj(y) / N

  S <- rep(0, n_freq)
  FF <- rep(0, n_freq)
  amp <- rep(NA, n_freq)
  A <- matrix(0, n_freq, n_freq)
  Q_matrix <- matrix(0, n_freq, n_freq)
  f <- list()
  B <- list()
  x <- list()
  freqs <- 2. * pi * seq(0, (N - 1)) / N
  x[[1]] <- x_data

  for (m in seq(n_freq)) {
    hx <- h_N * x[[m]]
    if (is.na(nu[m])) {
      f_base <- freqs[which.max(power(hx))]
      brackets <- c(f_base - pi / N, f_base + pi / N)
      brackets[1] <- max(min_freq, brackets[1])
      brackets[2] <- min(max_freq, brackets[2])

      tomax <- function(t) {
        function(f) {
          ft <- f * t
          Mod(hx %*% cis(-ft))
        }
      }
      f_max <- cmna::goldsectmax(tomax(t),
        brackets[1], brackets[2],
        tol = 1.e-10, m = 9999
      )
      if (abs(f_max) > abs(freqs[2] / 2)) {
         nu[m] <- f_max
      } else {
        nu[m] <- 0.
      }
    }
    f[[m]] <- cis(nu[m] * t)

    Q_matrix2 <- Q_matrix

    if (!fast) {
      for (i in seq(m)) {
        num <- (nu[m] - nu[i]) * N2
        Qm <- ifelse(num == 0, 1, Q(num))

        Q_matrix[m, i] <- cis(num) * Qm
        Q_matrix[i, m] <- Conj(Q_matrix[m,i])
      }
    } else {
      for (i in seq(m)) {
        Q_matrix[m, i] <- h_prod(f[[m]], f[[i]])
        Q_matrix[i, m] <- Conj(Q_matrix[m,i])
    }
  }

  A[m, ] <- 0
  if (m > 1) {
    f_m_bi <- rep(0, (m - 1))
    # eq. 17
    for (j in seq(m - 1)) for (s in seq(j)) f_m_bi[j] <- f_m_bi[j] - A[j, s] * (Q_matrix[m,s])       # eq. 19
  }

  norm <- 1
  # norm2 <- 1

  if (m>1) {
    norm <- norm + sum((f_m_bi)*Conj(f_m_bi))
    #for (j in seq(m-1)) norm <- norm -  2 * Re ( sum(A[j,(1:j)] * Q_matrix[j,m])  )
    for (j in seq(m-1)) norm <- norm -  2 * Re ( f_m_bi[j] *  sum(A[j,(1:j)] * Q_matrix[(1:j),m])  )
  }
  # if (m > 1) for (j in seq(1, (m - 1)))  norm2 <- norm2 - Mod(f_m_bi[j])^2
  # print (sprintf("norm1 = %.4f, norm2 = %.4f", norm, norm2))
  A[m,m ] <- 1. / sqrt(norm)

  if (m>1) for (j in seq(m - 1)) for (s in seq(j, (m - 1))) {
    A[m, j] <- A[m, j] + A[m,m] * ( f_m_bi[s] ) * (A[s, j]  )
  }

  # for test only 
  B[[m]] = 0. 
  # print (sprintf("A[%i,%i] = %.4f", m,m,A[m,m]))
  for (k in seq(m)) B[[m]] = B[[m]] + A[m,k] * f[[k]]   * cis(nu[k] * N2)
  # for (k in seq(m)) {
  #     print(sprintf("crossprod0 B[[%i]] * B[[%i]] = %.8f", k, m, Mod(h_prod(B[[k]],B[[m]]))))
  #     print(sprintf("crossprod2 B[[%i]] * B[[%i]] = %.8f", k, m, {
  #                                    crossprod <- 0;
  #                                    for (j in (seq(k))) for (jp in (seq(m))) crossprod <- crossprod + A[k,j] * Conj(A[m,jp]) * Q_matrix[j,jp] * cis((nu[j]-nu[jp])*N2);
  #                                    Mod(crossprod) }))
  # }
  #
  FF[m] <- h_prod(x[[m]], f[[m]])
  S[m] <- A[m,m] * FF[m] * cis(nu[m] * N2)
  

  #Stest = 0
  #for (j in seq(m)) Stest <- Stest + FF[j] * A[m, j] * cis(nu[j] * N2)
  # print (sprintf("S[m] = %.4f  + 1i %.4f and  Stest = %.4f + 1i * %.4f", Re(S[m]), Im(S[m]), Re(Stest), Im(Stest)))

  x[[m + 1]] <- x[[m]]
  for (j in seq(m)) x[[m + 1]] <- x[[m + 1]] - S[m] * A[m, j] * f[[j]] * cis(-(nu[j]) * N2)
  # these two lines are equivalent
  #for (j in seq(m)) x[[m + 1]] <- x[[m + 1]] - A[m,m] * FF[m] * A[m,j] * f[[j]] * cis((nu[m] - nu[j]) * N2)
  }

  m_max <- n_freq

  amp[1:m_max] <- 0
  for (s in seq(m_max)) for (j in seq(s,m_max)) amp[s] <- amp[s] + A[j,j]*A[j,s]*FF[j]*cis((nu[j]-nu[s]) * N2)

  phase <- Arg(amp)
  amp <- Mod(amp)
  OUT <- data.frame(nu = nu, amp = amp, phase = phase)
  return(OUT)
  return(OUT)
}




#' Frequency Modified Fourier transform for real and  complex series
#'
#' Iteratively extracts dominant frequencies from a real-valued time series
#' using a modified Fourier transform with frequency correction. The function
#' applies a Hanning window to the data, and amplitudes are renormalised through
#' a Gram–Schmidt procedure based on cross-products of Fourier components.
#' The cross-products can be computed either via analytical approximation
#' (faster, following Sidlichovsky and Nesvorny) or explicitly by numerical
#' integration.
#'
#' If frequencies \code{nu} are supplied, the routine uses them directly;
#' otherwise, it iteratively maximises the Fourier transform using golden
#' section search. The
#' function is typically called from [mfft_real()] once, and possibly a second
#' time to refine amplitudes and frequencies iteratively.
#'
#' @importFrom cmna goldsectmax
#' @param x_data The data provided either as a time series (advised), or as a vector.
#' @param min_freq,max_freq If provided, bracket the angular frequencies to be probed.
#'        Note: these are angular velocities (\eqn{2\pi / \mathrm{period}}), expressed in time-inverse units,
#'        with the time resolution encoded in `x_data` if the latter is a time series.
#'        The default for `min_freq` is 0, and for `max_freq` is \eqn{\pi}.
#' @param correction  0: no frequency correction (equivalent to Laskar);
#'        1: frequency correction using linear approximation;
#'        2: frequency correction using synthetic data;
#'        3: second order-correction using synthetic data (all documented in the Sidlichovsky and Nesvorny reference)
#' @param n_freq The number of frequencies returned, must be smaller than the length of `x_data`.
#' @param fast (default = TRUE) Uses analytical formulations for the crossproducts involving sines and cosines.
#'        Note: this is not really faster because the bottleneck is actually the golden section search,
#'        but more elegant.
#' @return a `discreteSpectrum` object, based on a data.frame with columns "Freq", "Amp", and "Phases".
#' @author Michel Crucifix
#' @references
#' \insertRef{sidlichovsky97aa}{gtseries}
#' @examples
#'
#' data(harmonic_sample)
#' data_to_analyse <- develop(harmonic_sample)
#' spec <- mfft_real(data_to_analyse)
#' print(spec)
#'
#' @rdname mfft
#' @export mfft_complex
mfft_complex <- function(x_data, n_freq = 10, min_freq = NULL, max_freq = NULL, correction = 1, fast = TRUE) {
  if (correction == 3) "this correction scheme is currently not implemented for real time series"
  N <- length(x_data)
  N2 <- N / 2.
  x_data <- stats::as.ts(x_data)
  dt <- deltat(x_data)

  my_min_freq <- ifelse(is.null(min_freq), 0, min_freq * dt)
  my_max_freq <- ifelse(is.null(max_freq), 2*pi, max_freq * dt)

  start_x <- stats::start(x_data)[1]
  N <- length(x_data)
  OUT <- mfft_complex_analyse(x_data, n_freq, fast, NULL, my_min_freq, my_max_freq)
 print(OUT)
  if (correction == 2) {
    x_data_synthetic <- rep(0, N)
    t <- seq(N) - 1
    for (i in seq(n_freq)) x_data_synthetic <- x_data_synthetic + OUT$amp[i] * cis(OUT$nu[i] * t + OUT$phase[i])
    OUT2 <- mfft_complex_analyse(x_data_synthetic, n_freq, fast, NULL, my_min_freq, my_max_freq)
    OUT$nu <- OUT$nu + (OUT$nu - OUT2$nu)
    OUT$amp <- OUT$amp + (OUT$amp - OUT2$amp)
    OUT$phase <- OUT$phase + (OUT$phase - OUT2$phase)
  } else if (correction == 1) {
    for (j in seq(n_freq)) {
      epsilon <- 0
      if ((j + 1) <= n_freq) {
        for (s in seq(j + 1, n_freq)) {
          print((OUT$nu[s] - OUT$nu[j]) * N2) * cos((OUT$nu[j] - OUT$nu[s]) * N2 + OUT$phase[j] - OUT$phase[s])
          epsilon <- epsilon + OUT$amp[s] *
            (
              Q_prime((OUT$nu[s] - OUT$nu[j]) * N2) * cos((OUT$nu[j] - OUT$nu[s]) * N2 + OUT$phase[j] - OUT$phase[s])
            )
        }
      }
      epsilon <- epsilon / Q_second_0 / N2 / OUT$amp[j]
      print(sprintf("correction = %.f5", epsilon))

      OUT$nu[j] <- OUT$nu[j] - epsilon
    }
    OUT <- mfft_complex_analyse(x_data, n_freq, fast, nu = OUT$nu, my_min_freq, my_max_freq)
  }

  # fold frequencies
  OUT$nu[which (OUT$nu > pi)] <- OUT$nu[which (OUT$nu > pi)] - (2*pi)
  
  OUT$nu <- OUT$nu / dt
  OUT$phase <- OUT$phase - start_x * OUT$nu

  to_be_corrected <- which(OUT$amp < 0)
  if (length(to_be_corrected)) {
    OUT$amp[to_be_corrected] <- -OUT$amp[to_be_corrected]
    OUT$phase[to_be_corrected] <- OUT$phase[to_be_corrected] + pi
  }

#  to_be_corrected <- which(OUT$nu < 0)
#  if (length(to_be_corrected)) {
#    OUT$nu[to_be_corrected] <- -OUT$nu[to_be_corrected]
#    OUT$phase[to_be_corrected] <- -OUT$phase[to_be_corrected]
#  }

  o <- order(OUT$amp, decreasing = TRUE)
  OUT$amp <- OUT$amp[o]
  OUT$nu <- OUT$nu[o]
  OUT$phase <- OUT$phase[o]

  OUT$phase <- (OUT$phase + (2 * pi)) %% (2 * pi)

  names(OUT) <- c("Freq", "Amp", "Phases")
  class(OUT) <- c("discreteSpectrum", "data.frame")
  attr(OUT, "data") <- x_data
  attr(OUT, "n_freq") <- n_freq
  return(OUT)
}

# Unified documentation ----------------------------------------------------

#' Modified Fourier Transform with frequency correction
#'
#' R implementation of the Modified Fourier Transform (MFT) based on the 
#' algorithm of Sidlichovský and Nesvorný (1997). The implementation prioritises
#' code clarity, portability, and extensibility over efficiency: it is slower
#' than the original C implementation, but straightforward to adapt to other
#' languages and contexts.
#'
#' Two variants are provided:
#' \itemize{
#'   \item \code{mfft_real}: designed for real-valued time series (original development).
#'   \item \code{mfft_complex}: designed for complex-valued time series.
#'   \item \code{mfft}: sends to one or the other depending on the data
#' }
#'
#' Both routines iteratively extract dominant frequencies, applying a Hanning
#' window and Gram–Schmidt orthogonalisation. In `mfft_complex`, the signal
#' is expressed as sum of complex exponentiasls; in `mfft_real`,it is expressed as a sum
#' of sinusoids (cosines and sines). Both routines allow for optional frequency correction at
#' different levels. The process starts with identifying the Fourier component
#' corresponding to the maximum power, storing the cosine of the corresponding
#' frequency as an element of the decomposition, and orthogonalising. In `mffht_real`,
#' a cosine function is taken first, followed by a corresponding sine function is then included, followed by further
#' orthogonalisation. Compared with the exponential basis \eqn{e^{i \omega t}},
#' the treatment for real signals attention to the non-unit norm of cosines and sines, and an
#' adapted frequency correction procedure.
#'
#' @param x_data Numeric vector or time series (real or complex, depending on
#'   the routine). The signal to analyse.
#' @param n_freq Integer. Number of frequencies to extract.
#' @param min_freq,max_freq Optional numeric scalars. Bounds for the angular
#'   frequencies considered in the search. Units depend on the sampling of
#'   \code{x_data}. Defaults are 0 and \eqn{\pi} for real series, and
#'   \eqn{-\pi} and \eqn{\pi} for complex series.
#' @param correction Frequency correction scheme (integer code). Options are:
#'   \itemize{
#'     \item 0: no correction (equivalent to Laskar).
#'     \item 1: linear approximation.
#'     \item 2: correction using synthetic data.
#'     \item 3: second-order correction (currently implemented for complex only).
#'   }
#' @param fast Logical. If \code{TRUE}, uses analytical formulations for the
#'   cross-products of Fourier components; if \code{FALSE}, they are computed
#'   numerically.
#'
#' @return A \code{discreteSpectrum} object, based on a \code{data.frame} with
#'   columns:
#'   \itemize{
#'     \item \code{Freq}: estimated angular frequencies
#'     \item \code{Amp}: amplitudes of the extracted components
#'     \item \code{Phases}: phases of the extracted components
#'   }
#'
#' @note The adaptation of the Modified Fourier Transform for real-valued time
#'   series (\code{mfft_real}) is an original development of the author.
#'
#' @references
#' Sidlichovský, M., & Nesvorný, D. (1997). Frequency modified Fourier transform
#' method. \emph{Celestial Mechanics and Dynamical Astronomy}, 65(2–3), 137–148.
#'
#' @author Michel Crucifix
#' @aliases mfft_real mfft_complex mfft
#' @rdname mfft
#' @export
mfft_real <- function(x_data, n_freq = 5, min_freq = NULL, max_freq = NULL, correction = 1, fast = TRUE) {
  if (correction == 3) "this correction scheme is currently not implemented for real time series"
  N <- length(x_data)
  N2 <- N / 2.
  x_data <- stats::as.ts(x_data)
  dt <- deltat(x_data)

  my_min_freq <- ifelse(is.null(min_freq), 0, min_freq * dt)
  my_max_freq <- ifelse(is.null(max_freq), pi, max_freq * dt)

  start_x <- stats::start(x_data)[1]
  N <- length(x_data)
  OUT <- mfft_real_analyse(x_data, n_freq, fast, NULL, my_min_freq, my_max_freq)

  if (correction == 2) {
    x_data_synthetic <- rep(0, N)
    t <- seq(N) - 1
    for (i in seq(n_freq)) x_data_synthetic <- x_data_synthetic + OUT$amp[i] * cos(OUT$nu[i] * t + OUT$phase[i])
    OUT2 <- mfft_real_analyse(x_data_synthetic, n_freq, fast, NULL, my_min_freq, my_max_freq)
    OUT$nu <- OUT$nu + (OUT$nu - OUT2$nu)
    OUT$amp <- OUT$amp + (OUT$amp - OUT2$amp)
    OUT$phase <- OUT$phase + (OUT$phase - OUT2$phase)
  } else if (correction == 1) {
    for (j in seq(n_freq)) {
      epsilon <- OUT$amp[j] * Q_prime(-2 * OUT$nu[j] * N2) * cos(2 * OUT$nu[j] * N2 + 2 * OUT$phase[j])
      if ((j + 1) <= n_freq) {
        for (s in seq(j + 1, n_freq)) {
          epsilon <- epsilon + OUT$amp[s] *
            (
              Q_prime((OUT$nu[s] - OUT$nu[j]) * N2) * cos((OUT$nu[j] - OUT$nu[s]) * N2 + OUT$phase[j] - OUT$phase[s]) -
                Q_prime((OUT$nu[s] + OUT$nu[j]) * N2) * cos((OUT$nu[j] + OUT$nu[s]) * N2 + OUT$phase[j] + OUT$phase[s])
            )
        }
      }
      epsilon <- epsilon / Q_second_0 / N2 / OUT$amp[j]

      OUT$nu[j] <- OUT$nu[j] - epsilon
    }
    OUT <- mfft_real_analyse(x_data, n_freq, fast, nu = OUT$nu, my_min_freq, my_max_freq)
  }

  OUT$nu <- OUT$nu / dt
  OUT$phase <- OUT$phase - start_x * OUT$nu

  to_be_corrected <- which(OUT$amp < 0)
  if (length(to_be_corrected)) {
    OUT$amp[to_be_corrected] <- -OUT$amp[to_be_corrected]
    OUT$phase[to_be_corrected] <- OUT$phase[to_be_corrected] + pi
  }

  to_be_corrected <- which(OUT$nu < 0)
  if (length(to_be_corrected)) {
    OUT$nu[to_be_corrected] <- -OUT$nu[to_be_corrected]
    OUT$phase[to_be_corrected] <- -OUT$phase[to_be_corrected]
  }

  o <- order(OUT$amp, decreasing = TRUE)
  OUT$amp <- OUT$amp[o]
  OUT$nu <- OUT$nu[o]
  OUT$phase <- OUT$phase[o]

  OUT$phase <- (OUT$phase + (2 * pi)) %% (2 * pi)

  names(OUT) <- c("Freq", "Amp", "Phases")
  class(OUT) <- c("discreteSpectrum", "data.frame")
  attr(OUT, "data") <- x_data
  attr(OUT, "n_freq") <- n_freq
  return(OUT)
}

#' @rdname mfft
#' @export
mfft_complex <- function(x_data, n_freq = 10, min_freq = NULL, max_freq = NULL, correction = 1, fast = TRUE) {
  if (correction == 3) "this correction scheme is currently not implemented for real time series"
  N <- length(x_data)
  N2 <- N / 2.
  x_data <- stats::as.ts(x_data)
  dt <- deltat(x_data)

  my_min_freq <- ifelse(is.null(min_freq), -pi, min_freq * dt)
  my_max_freq <- ifelse(is.null(max_freq), pi, max_freq * dt)

  start_x <- stats::start(x_data)[1]
  N <- length(x_data)
  OUT <- mfft_complex_analyse(x_data, n_freq, fast, NULL, my_min_freq, my_max_freq)

  if (correction == 2) {
    x_data_synthetic <- rep(0, N)
    t <- seq(N) - 1
    for (i in seq(n_freq)) x_data_synthetic <- x_data_synthetic + OUT$amp[i] * cis(OUT$nu[i] * t + OUT$phase[i])
    OUT2 <- mfft_complex_analyse(x_data_synthetic, n_freq, fast, NULL, my_min_freq, my_max_freq)
    OUT$nu <- OUT$nu + (OUT$nu - OUT2$nu)
    OUT$amp <- OUT$amp + (OUT$amp - OUT2$amp)
    OUT$phase <- OUT$phase + (OUT$phase - OUT2$phase)
  } else if (correction == 1) {
    for (j in seq(n_freq)) {
      epsilon <- 0
      if ((j + 1) <= n_freq) {
        for (s in seq(j + 1, n_freq)) {
          print((OUT$nu[s] - OUT$nu[j]) * N2) * cos((OUT$nu[j] - OUT$nu[s]) * N2 + OUT$phase[j] - OUT$phase[s])
          epsilon <- epsilon + OUT$amp[s] *
            (
              Q_prime((OUT$nu[s] - OUT$nu[j]) * N2) * cos((OUT$nu[j] - OUT$nu[s]) * N2 + OUT$phase[j] - OUT$phase[s])
            )
        }
      }
      epsilon <- epsilon / Q_second_0 / N2 / OUT$amp[j]
      print(sprintf("correction = %.f5", epsilon))

      OUT$nu[j] <- OUT$nu[j] - epsilon
    }
    OUT <- mfft_complex_analyse(x_data, n_freq, fast, nu = OUT$nu, my_min_freq, my_max_freq)
  }

  OUT$nu <- OUT$nu / dt
  OUT$phase <- OUT$phase - start_x * OUT$nu

  to_be_corrected <- which(OUT$amp < 0)
  if (length(to_be_corrected)) {
    OUT$amp[to_be_corrected] <- -OUT$amp[to_be_corrected]
    OUT$phase[to_be_corrected] <- OUT$phase[to_be_corrected] + pi
  }

#  to_be_corrected <- which(OUT$nu < 0)
#  if (length(to_be_corrected)) {
#    OUT$nu[to_be_corrected] <- -OUT$nu[to_be_corrected]
#    OUT$phase[to_be_corrected] <- -OUT$phase[to_be_corrected]
#  }

  o <- order(OUT$amp, decreasing = TRUE)
  OUT$amp <- OUT$amp[o]
  OUT$nu <- OUT$nu[o]
  OUT$phase <- OUT$phase[o]

  OUT$phase <- (OUT$phase + (2 * pi)) %% (2 * pi)

  names(OUT) <- c("Freq", "Amp", "Phases")
  class(OUT) <- c("discreteSpectrum", "data.frame")
  attr(OUT, "data") <- x_data
  attr(OUT, "n_freq") <- n_freq
  return(OUT)
}

#' @rdname mfft
#' @export mfft
#' @param force_complex in the `mfft` function, assume x_data is complex, even if it contains real numbers only
mfft <- function(x_data, n_freq = 15, min_freq = NULL, max_freq = NULL, correction = 1, force_complex = FALSE) {
  if (is.complex(x_data) || force_complex) {
    return(mfft_complex(x_data, n_freq, min_freq, max_freq, correction))
  } else {
    return(mfft_real(x_data, n_freq, min_freq, max_freq, correction))
  }
}


