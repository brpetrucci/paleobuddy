#' Bin occurrences
#'
#' \code{binner} takes the time interval vector in question and the list of
#' occurrences for a species.
#'
#' @param x the list containing occurrence times for a given species.
#'
#' @param bins a vector of time intervals corresponding to geological time ranges.
#'
#' @return a list of upper bounds for the occurrence times of a species.
#'
#' @author written by Matheus Januario.
#'
#' @examples
#'
#' # this function is simple to test: let us just create an artificial bin and
#' # occurrence lists and then check it bounds them correctly
#' x1 <- c(5.2, 4.9, 4.1, 3.2, 1, 0.2)
#' bins1 <- c(6, 5, 4, 3, 2, 1, 0)
#' binneredSamp1 <- binner(x1, bins1) # should be (1, 2, 1, 0, 1, 1)
#'
#' # it should work with any type of number in IntVec
#' x2 <- c(6.7, 5.03, 4.2, 3.4, 1.2, 0.4)
#' bins2 <- c(7, bins1)
#' binneredSamp2 <- binner(x2, bins2) # should be (1, 1, 1, 1, 0, 1, 1)
#'
#' # let us try with a real simulated species fossil record
#' sim <- BDSim(1, 0.1, 0.09, 15)
#' sampled1 <- Sample(1, sim$TE, sim$TS, 1, 15)
#' bins3 <- c(15.1, 12.3, 10, 7.1, 5.8, 3.4, 2.2, 0)
#' binneredSample <- binner(sampled1, bins3)
#' # one can then check by inspection that the binner is indeed counting correctly
#'
#' @name binner
#' @rdname binner
#' @export

binner<-function(x, bins){
  # create result
  res<-vector()

  # in case time was passed wrong
  bins <- sort(bins, decreasing=T)
  x <- sort(x, decreasing=T)

  # for each bin,
  for(i in 2:(length(bins)-1)){
    # get the occurrences before this bin and after previous
    res<-c(res, sum(x<bins[i-1] & x>=bins[i]))
  }

  # find the last one
  i<-length(bins)
  res<-c(res, sum(x<bins[i-1] & x>=bins[i]))
  return(res)
}
