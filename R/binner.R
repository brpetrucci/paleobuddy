#' Bin occurrences in geologic intervals
#'
#' Given a vector of fossil occurrences and time bins to represent geological
#' ranges, returns the occurrence counts in each bin.
#'
#' @param x The vector containing occurrence times for one given species.
#'
#' @param bins A vector of time intervals corresponding to geological time 
#' ranges.
#'
#' @return A vector of occurrence counts for each interval, sorted from 
#' furthest to closest to zero.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci
#'
#' @details The convention is, for each bin, to include all occurrences exactly
#' in the boundary furthest from zero and exclude bins exactly in the boundary
#' closest to zero. Then, in the bin closest to zero (i.e., the "last", or 
#' "most recent" bin), include all occurrence on each of the two boundaries.
#'
#' @examples
#'
#' ###
#' # first let us create some artificial occurrence data and check
#' 
#' # occurrence vector
#' x <- c(5.2, 4.9, 4.1, 3.2, 1, 0.2)
#' 
#' # bins vector
#' bins <- c(6, 5, 4, 3, 2, 1, 0)
#' 
#' # result
#' binnedSamp <- binner(x, bins)
#' binnedSamp
#' 
#' ###
#' # it should work with any type of number in bins
#' 
#' # occurrence vector
#' x <- c(6.7, 5.03, 4.2, 3.4, 1.2, 0.4)
#' 
#' # bins vector
#' bins <- c(7.2, 6.1, 5.6, 4.3, 3.2, sqrt(2), 1, 0)
#' 
#' # result
#' binnedSamp <- binner(x, bins)
#' binnedSamp
#' 
#' ###
#' # let us try with a real simulated species fossil record
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim(1, lambda = 0.1, mu = 0.05, tMax = 15)
#' 
#' # set seed
#' 
#' 
#' # sample it
#' sampled <- sample.species(sim = sim, rho = 1, tMax = 15, S = 1)
#' 
#' # bins vector
#' bins <- c(15.1, 12.3, 10, 7.1, 5.8, 3.4, 2.2, 0)
#' 
#' # result
#' binnedsample <- binner(sampled, bins)
#' binnedsample
#' 
#' @name binner
#' @rdname binner
#' @export

binner <- function(x, bins) {
  # check that bins and x are numeric vectors
  if (!is.numeric(x) || !is.numeric(bins)) {
    stop("x and bins must be numeric vectors")
  }
  
  # check there is more than 1 bin
  if (length(bins) < 2) {
    stop("bins requires at least two time points")
  }
  
  # create result
  res <- vector()
  
  # in case time was passed wrong
  bins <- sort(bins, decreasing = TRUE)
  x <- sort(x, decreasing = TRUE)
  
  # for each bin (except last)
  for (i in 1:(length(bins) - 1)) {
    # get the occurrences between this bin and the previous
    res <- c(res, sum(x <= bins[i] & x > bins[i + 1]))
  }
  
  # add bins exactly at the smaller bin boundary
  res[length(res)] <- res[length(res)] + sum(x == bins[length(bins)])
  
  return(res)
}

