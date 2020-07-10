#' Bin occurrences in geologic intervals
#'
#' \code{binner} takes the time interval vector in question and the list of
#' occurrences for a species.
#'
#' @param x the list containing occurrence times for a given species.
#'
#' @param bins a vector of time intervals corresponding to geological time ranges.
#'
#' @return a list of occurrence counts for each interval.
#'
#' @author written by Matheus Januario and Bruno do Rosario Petrucci
#'
#' @examples
#'
#' # this function is simple to test: let us just create an artificial bin and
#' # occurrence lists and then check it bounds them correctly
#' 
#' # occurrence list
#' x <- c(5.2, 4.9, 4.1, 3.2, 1, 0.2)
#' 
#' # bins list
#' bins <- c(6, 5, 4, 3, 2, 1, 0)
#' 
#' # result
#' binnedSamp <- binner(x, bins)
#' 
#' # check
#' print(sum(binnedSamp == c(1, 2, 1, 0, 1, 1)) == length(binnedSamp))
#' 
#' # it should work with any type of number in bins
#' 
#' # occurrence list
#' x <- c(6.7, 5.03, 4.2, 3.4, 1.2, 0.4)
#' 
#' # bins list
#' bins <- c(7, 6, 5, 4, 3, 2, 1, 0)
#' 
#' # result
#' binnedSamp <- binner(x, bins)
#' 
#' # check
#' print(sum(binnedSamp == c(1, 1, 1, 1, 0, 1, 1)) == length(binnedSamp))
#' 
#' # let us try with a real simulated species fossil record
#' 
#' # run the simulation
#' sim <- BDSim(1, pp = 0.1, qq = 0.09, tMax = 15)
#' 
#' # sample it
#' sampled <- Sample(S = 1, TE = sim$TE, TS = sim$TS, rr = 1, tMax = 15)
#' 
#' # bins list
#' bins <- c(15.1, 12.3, 10, 7.1, 5.8, 3.4, 2.2, 0)
#' 
#' # result
#' binnedSample <- binner(sampled, bins)
#' # one can then check by inspection that the binner is indeed counting correctly
#'
#' @name binner
#' @rdname binner
#' @export

binner <- function(x, bins) {
  # create result
  res <- vector()

  # in case time was passed wrong
  bins <- sort(bins, decreasing = TRUE)
  x <- sort(x, decreasing = TRUE)

  # for each bin,
  for (i in 2:(length(bins) - 1)) {
    # get the occurrences before this bin and after previous
    res <- c(res, sum(x < bins[i-1] & x >= bins[i]))
  }

  # find the last one
  i <- length(bins)
  res <- c(res, sum(x < bins[i-1] & x >= bins[i]))
  return(res)
}
