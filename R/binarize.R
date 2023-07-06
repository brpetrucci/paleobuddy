#' Bin true occurrences into geologic intervals
#'
#' Given the output of \code{sample.clade(..., returnTrue = FALSE)}, returns 
#' the occurrence counts in each bin (i.e., the same as 
#' \code{sample.clade(..., returnTrue = TRUE)}). This helps to trace 
#' perfect parallels between both output formats of \code{sample.clade}.
#'
#' @param fossils A \code{data.frame} exactly as returned 
#' by sample.clade(..., returnTrue = FALSE). See \code{?sample.clade)} 
#' for details.
#'
#' @param bins A vector of time intervals corresponding to geological time 
#' ranges.
#'
#' @return A \code{data.frame} exactly as returned by
#' \code{sample.clade(..., returnTrue = TRUE)}. See \code{?sample.clade)} for
#' details.
#'
#' @author Matheus Januario.
#'
#' @details This function helps a user bin "true occurrences" directly into
#'  binned occurrences, allowing for comparisons among "perfectly known" 
#'  fossil records and records that have a certain resolution (given by the 
#'  \code{bins} parameter).
#'
#' @examples
#'
#' ###
#' # set seed
#' set.seed(1) 
#' 
#' # run a birth-death simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.05, tMax = 50)
#' 
#' # choose bins
#' bins <- seq(0, 50, by = 1)
#' 
#' # generate "true" fossil occurrences
#' fossils_true <- sample.clade(sim, rho = 1, tMax = 50, returnTrue = TRUE)
#' 
#' # bin the true occurrences
#' fossils_binned <- binarize(fossils_true, bins)
#' 
#' # compare
#' fossils_true
#' fossils_binned
#' 
#' @name binarize
#' @rdname binarize
#' @export
#' 

binarize <- function(fossils, bins) {
  # sort bins in decreasing order
  bins <- sort(bins, decreasing = TRUE)

  # start results data frame
  res <- fossils[-ncol(fossils)]
  
  # initialize max and min columns
  res$MaxT <- NA
  res$MinT <- NA
  
  # iterate through all fossils
  for (i in 1:nrow(fossils)) {
    # find the lowest bin with higher value
    id <- max(which(bins >= fossils$SampT[i])) 
    
    # set max to that bin
    res$MaxT[i] <- bins[id]
    
    # and min to the one right after
    res$MinT[i] <- bins[id + 1]
  }
  return(res) 
}

