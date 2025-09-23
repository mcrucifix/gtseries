#' Produce an evenly-spaced time series by simple interpolation
#'
#' Uses standard R routines 'approx' or 'spline' to generate an evenly-spaced time series
#' based on data supplied on a time coordinate grid. 
#' @param data The data input, typically a dataframe with two columns
#' @param tcoord The name of the column name of the time coordinate. Defaults to "x"
#' @param dcoord The name of the column name of the data coordinate. Defaults to "y"
#' @param scale Time scaling (the coordinate "tcoord" will be divided by `scale`). 
#' @param thin Discards data that are sampled with a time space less than 0.5 times the time resolution in the output
#' @param n (defaults to 2048) is the number of data points in the output
#' @param spline If TRUE, the standard `spline` routine will be used, otherwise, `approx`
#' @param xlim If provided, brackets the time output
#' @return a time series R object
#' @author Michel Crucifix
#' @importFrom Rdpack reprompt
#' @importFrom stats approx spline ts
#' @export approx_ts
approx_ts <- function(data, tcoord="x", dcoord="y", scale=1, n=2048, thin=FALSE, spline=FALSE, xlim=NULL) {
  x <- data[[tcoord]]
  y <- data[[dcoord]]
  
  if (!is.null(xlim)) {
    T <- which(x >= xlim[1] & x <= xlim[2])
    x <- x[T]
    y <- y[T]
  }
  
  deltat <- diff(range(x)) / n
  
  if (thin) {
    thin_f <- 0.5 * deltat
    tmp <- x
    repeat {
      dd <- abs(diff(tmp))
      drop <- dd < thin_f
      if (any(drop)) {
        tmp <- tmp[-(1 + which(drop)[1])]
      } else {
        break
      }
    }
    y <- y[match(tmp, x)]
    x <- x[match(tmp, x)]
  }
  
  # Interpolation
  if (spline) {
    dummy <- stats::spline(x, y, n = n)
  } else {
    dummy <- stats::approx(x, y, n = n)
  }
  
  ts_out <- stats::ts(dummy$y, start = dummy$x[1] / scale, deltat = diff(dummy$x[1:2]) / scale)
  return(ts_out)
}
