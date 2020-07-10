#' Constant and non-constant rate species sampling
#'
#' \code{Sample} takes a species number, a vector of speciation and extinction 
#' times, a sampling rate and a maximum time for simulation and returns a list
#' of occurrence times for each species.
#'
#' @param S the species number to be sampled. Since \code{Sample} will be called
#' by a wrapper using \code{lapply}, it is through \code{S} that we apply this
#' function.
#'
#' @param TE a vector of extinction times, usually an output of \code{BDSim}.
#'
#' @param TS a vector of speciation times, usually an output of \code{BDSim}.
#'
#' @param rr a sampling rate function. Can be created by \code{MakeRate} for
#' simplicity, but can be any time-varying function.
#'
#' @param tMax the maximum simulation time, used by \code{rexp_var}.
#'
#' @return a list of occurrences for that species.
#'
#' @author written by Bruno do Rosario Petrucci and Matheus Januario.
#'
#' @examples
#'
#' # Note: all examples use just 1 lineage and very large preservation rates just
#' # to be clearer to the reader. Most of the time preservation will never be 
#' # this high
#' 
#' # let us start with a linear increase in preservation rate
#' 
#' # simulate a group
#' sim <- BDSim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation was short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) { 
#'   sim <- BDSim(n0 = 1, pp = -.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # preservation function
#' r<-function(t) {
#'   return(50 + t*30)
#' }
#' 
#' # time
#' t<-seq(0, 10, by = 0.1)
#' 
#' # visualizing from the past to the present
#' plot(x = t, y = rev(r(t)), main="Simulated preservation", type = "l", 
#'      col = "red", xlab = "Mya", ylab = "preservation rate", 
#'      xlim = c(10, sim$TE[1]))
#' 
#' \dontrun{
#' # sample
#' occs <- Sample(S = 1, TS = sim$TS[1], TE = sim$TE[1], rr = r, tMax = 10)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, sim$TE[1]), 
#'      xlab = "Mya")
#' lines(t, rev(r(t)))
#' }
#' 
#' # now let us try a step function
#' 
#' # simulate a group
#' sim <- BDSim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation was short lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) { 
#'   sim <- BDSim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we can create the sampling rate here from a few vectors
#' 
#' # rates
#' rlist <- c(50, 20, 80)
#' 
#' # rate shift times -  this could be c(10, 6, 2) 
#' # and would produce the same function
#' rShifts <- c(0, 4, 8)
#' 
#' # create the rate to visualize it
#' r <- MakeRate(rlist, tMax = 10, fShifts = rShifts)
#' 
#' # time
#' t <- seq(0, 10, by = 0.1)
#' 
#' # visualizing the plot from past to present
#' plot(x = t, y = rev(r(t)), main = "Simulated preservation", type = "l", 
#'      col = "red", xlab = "Mya", ylab = "preservation rate", 
#'      xlim = c(10, sim$TE[1]))
#' 
#' \dontrun{
#' # sample
#' occs <- Sample(S = 1, TS = sim$TS[1], TE = sim$TE[1], rr = r, tMax = 10)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, sim$TE[1]), 
#'      xlab = "Mya")
#' 
#' # frontiers of each regime
#' abline(v = 10 - rShifts, col = "red") 
#' }
#' 
#' # we can create a step function in a different way as well
#' 
#' # simulate a group
#' sim <- BDSim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation was short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- BDSim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # preservation function
#' r<-function(t) {
#'   ifelse(t < 4, 50,
#'          ifelse(t < 8, 20,
#'                 80)) 
#' }
#' # note how this function should be exactly the same as the previous one
#' 
#' # time
#' t <- seq(0, 10, by = 0.1)
#' 
#' # visualizing the plot from past to present
#' plot(x = t, y = rev(r(t)), main = "Simulated preservation", type = "l", 
#'      col = "red", xlab = "Mya", ylab = "preservation rate", 
#'      xlim = c(10, sim$TE[1]))
#' 
#' \dontrun{
#' # sample 
#' occs <- Sample(S = 1, TS = sim$TS[1], TE = sim$TE[1], rr = r, tMax = 10)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, sim$TE[1]), 
#'      xlab = "Mya")
#' abline(v = 10 - rShifts, col = "red")
#' }
#' 
#' # one could also use an environmental function to generate a sampling rate,
#' # see MakeRate()
#' 
#' @name Sample
#' @rdname Sample
#' @export

Sample<-function(S, TE, TS, rr, tMax) {
  # invert times since simulation goes from 0 to tMax
  TE <- tMax - TE
  TS <- tMax - TS

  # start when the species was born, end when it died
  
  # TS for initial species is -0.01, make it 0
  Now<-ifelse(TS[S] < 0, 0, TS[S])
  
  # TE for extant species is tMax + 0.01, make it tMax
  End<-ifelse(TE[S] > tMax, tMax, TE[S])

  # make rr a function if it isn't
  r <- rr
  rr <- ifelse(is.numeric(r), Vectorize(function(t) r), r)

  # initialize vector
  sampled <- c()

  # while we have not reached the time of extinction
  while (Now < End) {
    # take the waiting time for sampling, using rexp_var()
    WaitTimeR <- ifelse(rr(Now) > 0, rexp_var(1, rr, Now, tMax), Inf)

    # advance in time
    Now <- Now + WaitTimeR

    # if sampling comes after extinction, we don't include this occurrence
    if (Now > End) break

    # append to the vector
    sampled <- c(sampled, Now)
  }

  # finally, we invert time so that it goes from tMax to 0
  sampled <- tMax - sampled

  return(sampled)
}
