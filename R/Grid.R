#' Returns a grided sample of occurrences for one species
#'
#' \code{Grid} takes the time interval vector in question and the list of 
#' occurrences for a species.
#'
#' @param \code{samp} the list containing the occurrence times for a given 
#' species.
#'
#' @param \code{IntVec} a vector of time intervals corresponding to 
#' geological time ranges. If \code{returnTrue} is false, \code{SampleClade} 
#' returns the member of \code{IntVec} right after the true occurrence time. 
#' In this way, we simulate the granularity in real world fossil records. 
#' If \code{returnTrue} is true, this is ignored.
#' 
#' @return a list of upper bounds for the occurrence times of a species.
#'
#' @author written by Bruno do Rosario Petrucci.
#'
#' @examples
#' 
#' this function is simple to test: let us just create an artificial grid and
#' occurrence lists and then check it bounds them correctly
samp1 <- c(0.2, 1, 3.2, 4.1, 4.9, 5.2)
IntVec1 <- c(0, 1, 2, 3, 4, 5, 6)
gridedSamp1 <- Grid(samp1, IntVec1) # should be (1, 2, 4, 5, 5, 6)
#'
#' it should work with any type of number in IntVec
IntVec2 <- c(0.4, 1.2, 3.4, 4.2, 5.03, 6.7)
gridedSamp2 <- Grid(samp1, IntVec2) # should be (0.4, 1.2, 3.4, 4.2, 5.03, 6.7)
#' 
#' let us try with a real simulated species fossil record
sim <- BDSim(1, 0.1, 0.09, 50)
sampled1 <- Sample(1, sim$TE, sim$TS, 1, 50)
IntVec3 <- c(0, 2.2, 3.4, 5.8, 7.1, 10, 12.3, 15.1)
gridedSample <- Grid(sampled1, IntVec3) 
#' one can then check by inspection that the grid represents the upper bounds

Grid<-function(samp, IntVec){
  # for each occurrence, find the interval that is closest to the occurrence
  # as an upper bound
  Grided<-lapply(1:length(samp),
                 function(y){IntVec[IntVec-samp[[y]]==
                                      min(IntVec[which(IntVec-samp[[y]]>0)]
                                          -samp[[y]])]})
  return(unlist(Grided))
}
