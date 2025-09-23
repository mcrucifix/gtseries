#' Maximum entropy method 
#'
#' Maximum entropy spectral analysis (MESA) is the statistical estimation of the power spectrum of a stationary time series using the maximum entropy (ME) method. The resulting ME spectral estimator is computed by fitting an autoregressive model to the input time series.
#'
#' @param A input, a time-series object
#' @param method Estimation method, passed to the standard `ar` function, which is one of "yule-walker", "burg", "ols", "mle", "yw". Defaults to "burg".
#' @param order.max maximum order of the autoregressive model being fitted
#' @param ... Additional arguments passed to \code{ar}.
#' @param x A \code{memObject}, the result from \code{mem()} (for plotting).
#' @param period Logical. If \code{TRUE}, the x axis is shown as period instead of frequency.
#' @param xaxp Arguments passed to adjust the x axis.
#' @param yaxp Arguments passed to adjust the y axis.
#' @return  a memObject
#' @references
#' \insertRef{Burg75aa}{gtseries}
#' \insertRef{percival98}{gtseries}
#' \insertRef{Ghil02aa}{gtseries}
#' \insertRef{pardoiguzquiza06aa}{gtseries}
#' \insertRef{pardo-iguzquiza21aa}{gtseries}
#' @importFrom graphics mtext 
#' @importFrom sfsmisc axis.period axis.frequency
#' @importFrom stats ar
#' @importFrom stats coef
#' @export  mem
#' @export plot.memObject
mem <- function (A,method='burg',order.max=90,...) {
  deltat<-deltat(A)
  minf  <- 2./length(A)

  f <- 10^(seq(log10(minf),log10(0.5),length.out=400))

  U <- stats::ar(A,method=method,order.max=order.max,...)

  ##  v.maice multiplied by deltat to have power density (I still need to really undersant
  ##    why, but this works)
  arc <- sapply (U$ar,function (x) {x})
  mem <- data.frame(frequency=f/deltat,power= arspec(f,arc)*U$var.pred*(deltat),row.names=NULL)

  attr(mem,"var") <- U$v.maice*(deltat)
  attr(mem,"L") <- U$order
  attr(mem,"class") <- c("memObject","data.frame")

  mem
}


#' @rdname mem
plot.memObject <- function (x,period=FALSE,xaxp=NULL,yaxp=NULL,...)
{

  plot(x$frequency,x$power,
       frame=FALSE,axes=FALSE,log="xy",type="l",xlab="",ylab="Power density",...)
  local({
  if (is.null(yaxp)) {yaxp <- par("yaxp"); yaxp[3] <- 1; par(yaxp=yaxp)}
  if (is.null(xaxp)) {xaxp <- par("xaxp"); xaxp[3] <- 1; par(xaxp=xaxp)}
   eaxis(2,outer.at=FALSE)
  
  if (period) {axis.period()} else {axis.frequency()}
  })
}

