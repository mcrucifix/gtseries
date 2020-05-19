# SSA as in Ghil et al., 
# Advances spectral methods for climatic time series
# Rev Geoph., 40 (1), art. 3 (2002)
# doi 10.1029/2000RG000092

# Has been validated and show to produce EXACTLY the same results
# as the Spectra toolkin on 25.9.08
#----------------------------------------------------------------------   
# using astronomical forcing as sample data

ssa <- function(X=X,N=length(X),M=M)
{

  # construct  Toeplitz Matrix for reference signal (eq. 6)
  # according to Vautard and Ghil, 1989

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
   
  A <- matrix(0,N-M+1,10)



   for (k in 1:10) { for (t in 1:(N-M+1))
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

  PCA <- matrix(0,10,N);


  for (t in 1:N) {for (k in 1:10)
   {
   PCA[k,t] <- PCA[k,t] + crossprod(A[t-Lt[t]:Ut[t]+1,k],rhok[Lt[t]:Ut[t],k]);
  }}

  LA <- vector("double",10)
  for (k in 1:10) {
   PCA[k,] <- PCA[k,]/Mt
   LA[k]=crossprod(PCA[k,],X[1:N])
  }

  SSAObject<- list(lambda=lambda,LA=LA,A=A,rhok=rhok,cij=cij,PCA=PCA)
  attr(SSAObject,"class") <- "SSAObject"
  SSAObject
  }



plot.SSAObject <- function (SSAObject,...)
{
  yat <- outer(seq(1:10),10^(seq(-3,5)),"*")
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


