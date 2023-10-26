# SSA as in Ghil et al., 
# Advances spectral methods for climatic time series
# Rev Geoph., 40 (1), art. 3 (2002)
# doi 10.1029/2000RG000092

# Has been validated and show to produce EXACTLY the same results
# as the Spectra toolkin on 25.9.08
#----------------------------------------------------------------------   
# using astronomical forcing as sample data

#' Singular spectrum analysis 
#'
#' @description Singular spectrum analysis as referred and described in Ghil et al.
#' and in the  \url{http://research.atmos.ucla.edu/tcd/ssa/}
#' results have been validated and shown to produce exactly the same results
#' as the Spectra toolking on 25.9.08

#' @param X the input series (typically a numeric)
#' @param N if provided,  the length of the input series to beconsiderd 
#' @param M the length of the sliding window
#' @param Nrec: number of components to be output
#' @importFrom sfsmisc axTexpr
#' @importFrom sfsmisc eaxis

#' @references
#' \insertRef{Ghil02aa}{gtseries}
#' \insertRef{VAUTARD89aa}{gtseries}
#' @return a SSA object, compatible with the supplied plot function
#' @examples
#' 
#' bbridge <- function (n=256, amp) {
#'   rw1 <- cumsum(rnorm(n,1))*amp / (n)
#'   rw2 <- cumsum(rnorm(n,1))*amp / (n)
#'   rw2 <- rw2 * rw1[n] / rw2[n]
#'   bbridge <- rw1 - rw2
#' }
#' 
#' N <- 256
#' 
#' Phase1 <- bbridge(N,2*pi)
#' Phase2 <- bbridge(N,2*pi)
#' Amp1 <- 1 + bbridge(N,0.5)
#' Amp2 <- 1 + bbridge(N,0.5)
#' 
#' t = seq(N)
#' 
#' P1 = 40
#' P2 = 90
#' 
#' sig1 <- sin(2*pi*t/P1 + Phase1)*Amp1
#' sig2 <- sin(2*pi*t/P2 + Phase2)*Amp2
#' 
#' noise <- bbridge(N, 0.3) + rnorm(N,0.2)
#' 
#' signal <- noise + sig1 + sig2
#' 
#' plot (signal, type='l')
#' 
#' M <- 50
#' 
#' SSAsignal <- ssa (signal, M=M, Nrec=M)
#' 
#' NRec <- c(10,30,M)
#' 
#' for (i in seq(along=NRec)) {
#' rec <- apply(SSAsignal$PCA[seq(NRec[i]),], 2, sum)
#' lines(rec,col=i+1)
#' }
#' 
#' # full reconstruction (modulo the mean of the signal)
#' fullrec <- apply(SSAsignal$PCA[seq(M),], 2, sum)
#' 
#' # show that the difference is small
#' plot(fullrec - signal + mean(signal))
#' 
#' @export ssa
ssa <- function(X=X,N=length(X),M=M,Nrec=10) {

  # construct  Toeplitz Matrix for reference signal (eq. 6)
  # according to Vautard and Ghil, 1989
  # added: Nrec : number of components of reconstruction

  X <- X - mean(X)
  cij <- vector("double",M)

  for (i in 1:M) {for (j in i:M)
     ij <- j-i
     T <- N-ij
     cij[ij+1]=crossprod(X[1:T],X[(1+ij):(T+ij)])/T
     }

  c <- toeplitz(cij)

  # calculate eigenvalues of reference signal (eq. 7)

  Eig <- eigen(c)
  lambda <- Eig$values  ## they are immediately sorted. Thanks R!
  rhok   <- Eig$vectors

  # sort Lambda of reference signal 

  #  [lambda,index]=sort(Lambdak);


  # evaluation of the principal components of the original signal (eq. 10)
   
  A <- matrix(0,N-M+1,Nrec)



   for (k in 1:Nrec) { for (t in 1:(N-M+1))
     A[t,k] <- crossprod(X[t:(t+M-1)],rhok[,k])
   }

  # reconstruction from the eigenvectors (eq. 11)

  Mt <- matrix(0,N,1)
  Lt <- matrix(0,N,1)
  Lt1 <- matrix(0,N,1)
  Ut <- matrix(0,N,1)
  Ut1 <- matrix(0,N,1)

  for (t in 1:(M-1)) {
    Mt[t]<-t;
    Lt[t]<-1;
    Ut[t]<-t}

  Mt[M:(N-M)]<-M
  Lt[M:(N-M)]<-1
  Ut[M:(N-M)]<-M

  for (t in (N-M+1):N) {
    Mt[t]<-N-t+1
    Lt[t]<-t-N+M
    Ut[t]<-M}

  PCA <- matrix(0,Nrec,N);


  for (t in 1:N) {for (k in 1:Nrec)
   {
   PCA[k,t] <- PCA[k,t] + crossprod(A[t-Lt[t]:Ut[t]+1,k],rhok[Lt[t]:Ut[t],k]);
  }}

  LA <- vector("double",Nrec)
  for (k in 1:Nrec) {
   PCA[k,] <- PCA[k,]/Mt
   LA[k]=crossprod(PCA[k,],X[1:N])
  }

  SSAObject<- list(lambda=lambda,LA=LA,A=A,rhok=rhok,cij=cij,PCA=PCA)
  attr(SSAObject,"class") <- "SSAObject"
  SSAObject
  }



#' @rdname ssa
#' @export
plot.SSAObject <- function (SSAObject,...)
{
  Nrec = length(SSAAbject$LA)
  yat <- outer(seq(1:Nrec),10^(seq(-3,5)),"*")
  ylab <- 10^(-2:5)

   ss <- lapply(seq(along = ylab),
                    function(i) if(ylab[i] == 0) quote(0) else
                       substitute(10^E, list(E=log10(ylab[i]))))

   S <- do.call("expression",ss)
  plot(SSAObject$lambda,frame=T,axes=F,log="y",xlab="rank",ylab="Eigenvalue",...)
  axis(1,xlab="rank")
  axis(2,at=yat,ylab="Eigenvalue",labels=F)
  axis(2,at=ylab,labels=axTexpr(1,at=ylab),tick=F,las=1)
}


