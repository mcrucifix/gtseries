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

    dummy <-  eval(call(ifelse(spline,"spline","approx"),x,y,n=n))
    out <<- ts(dummy$y,start=dummy$x[1]/scale, deltat=diff(dummy$x[1:2])/scale)
   })

  out
  }
