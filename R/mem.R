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


axis.period    <- function(side=1,...) { axis.log10(side,"Period",factor=-1,...) }
axis.frequency <- function(side=1,...) { axis.log10(side,"Frequency",factor=1,...) }


axis.log10 <- function(side=1,title="Power",factor=1,line=3,outer.at=TRUE,small.tcl=3/5,las=1,...)
{
   range10 <- sort(factor*par("usr")[c(3,4)-2*(side %% 2)])
   range <- 10^range10

   big_ticks <- 10^seq(floor(range10[1]),ceiling(range10[2]))
   big_ticks <- big_ticks[which(big_ticks > range[1] & big_ticks < range[2])]
   axis(side,at=big_ticks^factor,labels=axTexpr(side=side,at=big_ticks,drop.1=TRUE),las=las,...)

   if (outer.at ) {
         small_ticks <- outer(seq(1,9),10^seq(floor(range10[1]),ceiling(range10[2])),"*")
         small_ticks <- small_ticks[which(small_ticks > range[1] & small_ticks < range[2])]
         axis(side,at=small_ticks^factor,labels=FALSE,tcl=3/5*par("tcl"),...)
         }
   graphics::mtext(title,side=side,line=3,cex=par("cex.lab")*par("cex"))
}

add.intercept <- function (memObject, xlim = range(memObject$frequency)[2]*c(5.e-2,0.5),annote=TRUE)
{
frequency <- memObject$frequency
print(range(memObject$frequency)[2])

t <- which (frequency > xlim[1] & frequency < xlim[2])

fm1 <- lm(log10(power) ~ log10(frequency),as.data.frame(memObject[t,]))
lines(xlim,10^(stats::coef(fm1)[1]+log10(xlim)*stats::coef(fm1)[2]),lty=2)
#print(xlim)
#print(10^(coef(fm1)[1]+log10(xlim)*coef(fm1)[2]))
if (annote) {
text(mean(frequency)+0.1*diff(range(frequency)),mean(memObject$power)+0.0*diff(range(memObject$power)),
    paste("slope = ",format(coef(fm1)[[2]],digit=2)))}
coef(fm1)[[2]]
}
