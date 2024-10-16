#' Develop a spectrum into a time series (generic)
#' @param  arg: input class
#' @param  times: if supplied, times of the decomposition
#' @param  start: if supplied, overrides time and will generate a time series with start and deltat, which must then
#'         be supplied as well
#' @param  deltat : see start. 
#' @note   place holder for type-specific develop functions
#' @export develop
#' @return nothing
develop <- function(M, times=NULL, start=NULL, end=NULL, deltat=NULL,...){
     UseMethod("develop") 
}

