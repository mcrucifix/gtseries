#' MFFT reconstruction
#' @importFrom stats deltat start
#' @rdname mfft_deco
#' @param  M : mfft_deco object 
#' @param  sum : TRUE if user wants to sum components in the reconstruction
#' @export reconstruct_mfft
#' @return list of reconstructed components if sum=FALSE,  full
#'         reconstructed time series otherwise
reconstruct_mfft  <- function(M, sum=TRUE){
 if (!(attr(M,"class")=="mfft_deco")) stop ("object is not a MFFT decomposition")
 xdata <- attr(M,"data")
 nfreq <- attr(M,"nfreq")
 times <- seq(length(xdata))*stats::deltat(xdata) + stats::start(xdata)[1]
 reconstructed <- lapply(seq(nfreq), function(i) ts( M$Amp[i]*cos(M$Freq[i]*times + M$Phase[i]), start=stats::start(xdata), deltat=stats::deltat(xdata)) )

 if ( sum ) reconstructed <- ts(apply(simplify2array(reconstructed), 1 , sum), start=stats::start(xdata), deltat=stats::deltat(xdata))
}

#' MFFT ANOVA
#' not ready. do not use. 
#' @rdname mfft_deco
#' @export mfft_anova
mfft_anova  <- function(M){
 if (!(attr(M,"class")=="mfft_deco")) stop ("object is not a MFFT decomposition")
 xdata <- attr(M,"xdata")
 nfreq <- attr(M,"nfreq")
 N <- length(xdata)
 times <- seq(length(xdata))*deltat(xdata) + start(xdata)
 reconstructed <- sapply(seq(nfreq), function(i) M$Amp[i]*cos(M$Freq[i]*times + M$Phase[i]) )
 cum_reconstruct <- apply(reconstructed, 1, cumsum)
 residual_vars <- apply(apply(cum_reconstruct, 1, function(i) xdata-i) , 2, function(j) sum(j^2))
 var0 <- sum(xdata^2)
 master_vars <- c(var0, residual_vars[-length(residual_vars)])
 p2s <- 2*seq(nfreq)
 p1s <- c(0, p2s[-length(p2s)])
 F <- (master_vars - residual_vars)/residual_vars * (N - p2s)/(p2s-p1s)
}


#' @rdname mfft_deco
#' @export
as.data.frame.mfft_deco <- function(x) {data.frame(Freq=x$Freq, Amp=x$Amp, Phases=x$Phases)}


#' @rdname mfft_deco
#' @param a `mfft_deco` object, typically the output of a `mfft` call. 
#' @param labels to be set above the frequency peaks. Can be the output of `attributeTone`
#' @param periods if TRUE will add a lower axis with period labels
#' @export
plot.mfft_deco <- function (M,periods=FALSE,labels=NULL,...){
#   O <- order(M$Freq)
  plot(abs(M$Freq), abs(M$Amp),'h',ylab="Amplitudes", xlab="",  ...)
  if (periods) {
    frequencies <- pretty(range(M$Freq/(2*pi)))
    plabels <- as.character(1/frequencies)
    if (0 %in% frequencies) plabels[which(frequencies == 0)] = "∞"
    axis(1, line=3, at=2*pi*frequencies, labels=plabels)
    mtext("Rate", 1, 2)
    mtext("Period", 1, 4)
  } else {
    mtext("Rate", 1, 3)
  }
  points(abs(M$Freq), abs(M$Amp),'p',...)
  if (!is.null(labels)) {
    yshift <- 0.05*diff(range(M$Amp))
    text(M$Freq, M$Amp, labels, srt=90, , adj=-0.4)
  }
}

#' @rdname mfft_deco
#' @export
lines.mfft_deco <- function (M,...){
#   O <- order(M$Freq)
  lines(abs(M$Freq), abs(M$Amp),'h',...)
  points(abs(M$Freq), abs(M$Amp),'p',...)
}


#' @rdname mfft_deco
#' @export
print.mfft_deco <- function (M,...){
  print.data.frame(cbind(as.data.frame(M), Period=2*pi/M$Freq))
}



