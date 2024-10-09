#' Generation of combination of tones
#' 
#' Generates a vector with combinations of an input vector of frequencies, wih
#' explicit label names, up to order 3 (this could be made more flexible is the future)
#'
#' @importFrom RcppAlgos comboGeneral
#' @param omegas: vector of references frequencies, optionally with rownames, 
#' @param keepPositives : if TRUE, then only keeps positive combinations of frequencies
#' @return a vector with combination of tones and explicit rownames, using, if available, the
#'         rownames provided in the input vector omega
#' @author Michel Crucifix
#' @export toneCombinations
#' @examples
#' omegas <- c( 0.123, 0.14312, 0.33251, 0.554313)
#' print(toneCombinations(omegas))

toneCombinations <- function(omegas, keepPositives=TRUE){
 twoomegas <- c(-omegas,omegas)
 indices  <- c(-seq(length(omegas)), seq(length(omegas)))
 result = rbind(
             c(0, 0, 0), 
             cbind (   RcppAlgos::comboGeneral(twoomegas, m = 1, repetition = T), 0, 0), 
             cbind (   RcppAlgos::comboGeneral(twoomegas, m = 2, repetition = T) , 0), 
                       RcppAlgos::comboGeneral(twoomegas, m = 3, repetition = T))
 combos = rbind(
             c(0, 0, 0), 
             cbind (   RcppAlgos::comboGeneral(indices, m = 1, repetition = T), 0, 0), 
             cbind (   RcppAlgos::comboGeneral(indices, m = 2, repetition = T) , 0), 
                       RcppAlgos::comboGeneral(indices, m = 3, repetition = T))


 sumre  <- apply(result, 1, sum)
 whichUnique <- seq(length(sumre))[!duplicated(sumre)] 
 uniqueCombos <- combos[whichUnique, ]
 uniqueSums   <- sumre[whichUnique ]
 
 # further filtering (?)
 if (keepPositives){
   tokeep <- which(uniqueSums > 0)
   uniqueCombos <- uniqueCombos[tokeep,]
   uniqueSums <- uniqueSums[tokeep]
 }

 names(uniqueSums) <- apply(uniqueCombos, 1, function(g) generate_name(g,"s", labels=names(omegas)))
 return(cbind(uniqueSums))

}

generate_name <- function(invec,char="s", labels = NULL){
  # to do: if labels are supplied, to not produce
  tt <- as.data.frame(table(invec))
  flag <- FALSE
  tmp <- sapply(seq(nrow(tt)), function(i) {
    indice <- as.integer(as.character(tt[i,][[1]]))
    if (indice == 0) return("") else {
    if (flag) signchar <- ifelse(sign(indice)+1, " + "," - ") else 
               signchar <- ifelse(sign(indice)+1, "","-") 
    freq <- tt[i,][[2]]
    freqchar <-ifelse(freq == 1, "",  sprintf("%d",freq))
    flag <<- TRUE
    if (is.null(labels)) labelname <- sprintf("%s%d", char, abs(indice)) else
                         labelname <- labels[abs(indice)]
    return(sprintf("%s%s%s",signchar,freqchar,labelname))
                }})
#     Reduce(function(i) paste (i, sep=""), tmp)
    Reduce(function(i, j) paste (i, j, sep=""), tmp)
}


#' Attribution of combination of tones
#' 
#' Based on a vector of frequencies (`infreq`), and a vector of referenc
#' frequencies with row names (it will be input to `toneCombinations`), 
#' attribute the `infreq` frequencies with two possible degrees of tolerance
#' 
#' @param infreq : input frequencies
#' @param omegas : reference frequencies (a numeric vector which may contain explicit row names)
#' @param tol1 : acceptable tolerance for being considered as a certain attribution
#'               (if several frequencies match the criteria, the closest will be taken)
#' @param tol2 : acceptable tolerance for being considered as a likely or plausible
#' @export attributeTones
#'
#' @examples
#' omegas <- c( 0.123, 0.14312, 0.33251, 0.554313)
#' names(omegas) <- c('g1','g2','s1','s2')
#' outamps <- c(1., 2, 0.2 , 0.5, 0.5)
#' outfreqs <- c(1., 1.2432, omegas[1]+omegas[3]+0.00000002, omegas[1]-omegas[4]+0.00004, 0.15)
#' 
#' attributions <- attributeTones(outfreqs, omegas)
#' 
#' cbind(outfreqs, attributions)
#' 
#' plot(outfreqs, outamps, type='h')
#' text(outfreqs, outamps+0.1, attributions)
#'
attributeTones <- function(infreq , omegas, tol1 = 1.e-6, tol2 = 1.e-4) { 
  attributions <- rep("", length(infreq))
  combis <- toneCombinations(omegas)
  for (i in seq(infreq)){
    deltas <- abs(infreq[i] - combis)
    bestSuspect <- which.min(abs(infreq[i] - combis))
    if ( deltas[bestSuspect] < tol1 ) attributions[i] <- rownames(combis)[bestSuspect] else
      if ( deltas[bestSuspect] < tol2 ) attributions[i] <- paste(rownames(combis)[bestSuspect],"?")
  }
  return (attributions)
}

