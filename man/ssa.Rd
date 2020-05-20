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
