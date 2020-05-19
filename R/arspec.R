## Calculates the theoretical density frequency of an
## autoregressive process, according to equation
## (3.3) of Akaike (1969)

arspec <- function (f,ar){
Exp <- exp(-1i*2*pi*outer(f,seq_along(ar),"*"))
1./abs(1-Exp%*%ar)^2
}

