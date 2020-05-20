\name{ssa}
\alias{ssa}
\title{Singar Spectrum Analysis}
\description{
  Singular spectrum analysis based on Ghil et al. 2002}

\usage{
 (X = X, N = length(X), M = M, Nrec = 10)
}
\arguments{
  \item{X}{numeric vector}
  \item{N}{number of points to be considered in the vector}
  \item{M}{window length}
  \item{Nrec}{number of components to keep for the reconstruction}
}
\details{}
\value{
  an \code{SSAObject} 
}
\author{Michel Crucifix}
\references{
Yiou P., E. Baert and M. F. Loutre (1996), Spectral analysis of climate data, Surveys in Geophysics, (17) 619-663 doi:10.1007/BF01931784

Ghil M., M. R. Allen, M. D. Dettinger, K. Ide, D. Kondrashov, M. E. Mann, A. W. Robertson, A. Saunders, Y. Tian, F. Varadi and P. Yiou (2002), Advanced spectral methods for climatic time series, Reviews of Geophysics, (40) 769-794 doi:10.1029/2000RG000092
}

\examples{
require (gtseries)

bbridge <- function (n=256, amp) {
  rw1 <- cumsum(rnorm(n,1))*amp / (n)
  rw2 <- cumsum(rnorm(n,1))*amp / (n)
  rw2 <- rw2 * rw1[n] / rw2[n]
  bbridge <- rw1 - rw2
}

N <- 256

Phase1 <- bbridge(N,2*pi)
Phase2 <- bbridge(N,2*pi)
Amp1 <- 1 + bbridge(N,0.5)
Amp2 <- 1 + bbridge(N,0.5)

t = seq(N)

P1 = 40
P2 = 90

sig1 <- sin(2*pi*t/P1 + Phase1)*Amp1
sig2 <- sin(2*pi*t/P2 + Phase2)*Amp2

noise <- bbridge(N, 0.3) + rnorm(N,0.2)

signal <- noise + sig1 + sig2

plot (signal, type='l')

M <- 50

SSAsignal <- ssa (signal, M=M, Nrec=M)

NRec <- c(10,30,M)

for (i in seq(along=NRec)) {
rec <- apply(SSAsignal$PCA[seq(NRec[i]),], 2, sum)
lines(rec,col=i+1)
}

# full reconstruction (modulo the mean of the signal)
fullrec <- apply(SSAsignal$PCA[seq(M),], 2, sum)

# show that the difference is smal
plot(fullrec - signal + mean(signal))
}
