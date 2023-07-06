#' Bin true occurrences into geologic intervals
#'
#' Given the output of \code{sample.clade(..., returnTrue = FALSE)}, returns 
#' the occurrence counts in each bin (i.e., the same as 
#' \code{sample.clade(..., returnTrue = TRUE)}). This helps to trace 
#' perfect parallels between both output formats of \code{sample.clade()}.
#'
#' @param fossil A \code{data.frame} exactly as returned 
#' by sample.clade(..., returnTrue = FALSE). See \code{help(sample.clade)} 
#' for details.
#'
#' @param bins A vector of time intervals corresponding to geological time 
#' ranges.
#'
#' @return A \code{data.frame} exactly as returned 
#' by sample.clade(..., returnTrue = TRUE). See \code{help(sample.clade)} 
#' for details.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci
#'
#' @details This function helps a user to bins "true occurrences" directly into
#'  binned occurrences, allowing for comparisons among "perfectly known" fossils
#'  records and records that have a certain resolution (given by the \code{bins}
#'    inputted).
#'
#' @examples
#'
#' # run a birth-death simulation:
#' sim= bd.sim(n0 = 1, lambda = .1, mu = .05, tMax = 50)
#' 
#' # choose bins
#' bins = seq(0,50,by=1)
#' 
#' # generate "true" fossil occurrences
#' fsls1 = sample.clade(sim, rho = 1, tMax=50, returnTrue = T)
#' 
#' # bin the true occurrences:
#' fsls2 = binarize(fsls1, bins)
#' 
#' @name binarize
#' @rdname binarize
#' @export

binarize=function(fossil, bins){

  bins = sort(bins, decreasing = T)

  res=fossil[-ncol(fossil)]
  res$MaxT=NA
  res$MinT=NA
  for(i in 1:nrow(fossil)){
    id = hist(fossil$SampT[i], breaks=bins, plot=FALSE)$counts
    res$MaxT[i]= bins[which(id==1)+1]
    res$MinT[i]= bins[which(id==1)]
  }
  return(res) 
}

