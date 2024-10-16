#' Moving Modified Fourier Transform
#'
#' @param xdata The data provided either as a time series (advised), or as a vector. 
#' may be complex
#' @param seglength Length of moving segment. Both xdata length and seglenthm are adviced to be a power of 2 (faster)
#' @param ...  passed to mmft
#' @return a `mmfft` object
#' @note  in the current implementation, the right-hand-side of xdata that is left after an integer division 
#'        by seglenth is discarded
#' @author Michel Crucifix
#' @export mmfft
mmfft <- function(xdata, seglength = length(xdata) %/% 16, ...){
  N <- length(xdata)
  nsections <- (N %/% seglength) 
  N0 <- nsections * seglength
#   remainder  <- seglength - N0
  xdata0 <- matrix(  xdata[1:N0], byrow=FALSE, seglength, nsections)
  myMfft <- function(x) mfft(x, ...)
  OUT <- apply(xdata0, 2, myMfft)
  class(OUT) <- "mmfft"
  attr(OUT, "nsections") <- nsections
  attr(OUT, "seglength") <- seglength
  attr(OUT, "start") <- ifelse(is.ts(xdata), stats::start(xdata), 1)
  attr(OUT, "deltat") <- ifelse(is.ts(xdata), stats::deltat(xdata), 1)
  return(OUT)
}
  
#' @rdname mmfft
#' @export
plot.mmfft <- function(x){
  freqrange <- 
    c(min(sapply(x, function(xs) min(xs$Freq))), 
      max(sapply(x, function(xs) max(xs$Freq))))

  amprange <- 
    c(min(sapply(x, function(xs) min(xs$Amp))), 
      max(sapply(x, function(xs) max(xs$Amp))))

  amp2lwd <- function(amp){ 3*amp/amprange[2] }

  length <- attr(x, "nsections") * 
            attr(x, "seglength")

  nsec <- length(x)
  dt <- attr(x,"deltat")
  tsec <- attr(x,"seglength")*dt

  tstart <- attr(x,"start")
  tend   <- tstart  + nsec * tsec

  plot(c(tstart, tend), freqrange, type='n', xlab='Time', ylab='Rate')
 
  for (iseq in seq(nsec)){
    trange <- tstart + c(iseq-1,iseq)*tsec
    obj <- x[[iseq]]
    nfreq <- length(obj$Freq)
    lwds <- sapply(obj$Amp, amp2lwd)
    for (j in seq(nfreq)){
      lines(trange, rep(obj$Freq[j],2), lwd=lwds[j])
    }
  }

}

