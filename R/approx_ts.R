#' Produce an evenly-spaced time series by simple interpolation
#'
#' Uses standard R routine 'approx' or 'spline' to generate an evenly-spaced time series
#' based on data supplied on a time coordinate grid. 
#' @param data The data input, typically a  dataframe,
#'        with two columns
#' @param tcoord The name of the column name of the time coordinate. Defaults to "x"
#' @param dcoord The name of the column name of the data coordinate. Defaults to "y"
#' @param scale Time scaling (the coordinate "tcoord" will be divided by `scale`. 
#' @param thin Discards data that are sampled with a time space less than 0.5 times the time resolution in the output
#' @param  n (defaults to 2048) is the number of data points in the output

#' @param  spline : if True, the standard `spline` routine will be used, otherwise, `approx'
#' @param  xlim : if provided, brackets the time output
#' @return a time series R object
#' @author Michel Crucifix


#' @importFrom Rdpack reprompt
#' @export approx_ts
approx_ts <- function (data,tcoord="x",dcoord="y",scale=1,n=2048,thin=FALSE,spline=FALSE,xlim=NULL)
  {
  local({
   x <- data[tcoord][[1]] 
   y <- data[dcoord][[1]] 
   if (!is.null(xlim)) {T <- which(x>=xlim[1] & x <= xlim[2]) ; x<-x[T]; y <- y[T]}
   deltat <- diff(range(x))/n
   ## thin data according with thin factor deltat

   if (thin) {
   thin_f <- 0.5*deltat
   tmp <-x
   while (TRUE) { dd <- abs(diff(tmp))
    if (any(drop <- dd < thin_f))
      tmp <- tmp[-(1 + which(drop)[1])]
    else
      break
      }

    y <- y[match(tmp,x)]
    x <- x[match(tmp,x)]
    }

    ## then interpolates
    if (spline) {
      dummy <- spline(x, y, n = n)
    } else {
      dummy <- approx(x, y, n = n)
    }
    out <<- ts(dummy$y,start=dummy$x[1]/scale, deltat=diff(dummy$x[1:2])/scale)
   })

  out
  }
