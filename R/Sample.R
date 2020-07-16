#' Constant and time-dependent rate species sampling
#' 
#' Generates a list of occurrence times for a species in a simulation using a
#' constant or a function of absolute time as a rate for a Poisson process. For
#' sampling of more than one species and/or taking into account species age
#' instead of absolute time, see \code{sample.clade} and \code{sample.adpp}.
#' Note that while the Poisson process occurs in forward time, we return (both in
#' birth-death functions and here) results in backwards time, so that time is
#' inverted using \code{tMax} both at the beginning and end of \code{sample}.
#'
#' @param S The species number to be sampled. Since \code{sample} will be called
#' by a wrapper using \code{lapply}, it is through \code{S} that we apply this
#' function.
#'
#' @param sim A \code{sim} object, usually an output of \code{bd.sim}.
#'
#' @param rr A sampling rate function. Can be created by \code{make.rate} for
#' simplicity, but can be any time-varying function, or a constant.
#'
#' @param tMax The maximum simulation time, used by \code{rexp.var}. A sampling
#' time greater than \code{tMax} would mean the occurrence is sampled after the
#' present, so for consistency we require this argument. This is also required
#' to ensure time follows the correct direction both in the Poisson process and
#' in the return.
#'
#' @return A list of occurrences for that species.
#'
#' @author Bruno do Rosario Petrucci and Matheus Januario.
#'
#' @examples
#'
#' ###
#' # let us start with a linear increase in preservation rate
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation was short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # preservation function
#' r <- function(t) {
#'   return(1 + 0.25*t)
#' }
#' 
#' # time
#' t <- seq(0, 10, by = 0.1)
#' 
#' # visualizing from the past to the present
#' plot(x = t, y = rev(r(t)), main="Simulated preservation", type = "l",
#'      xlab = "Mya", ylab = "preservation rate",
#'      xlim = c(10, sim$TE[1]))
#' 
#' # sample
#' occs <- sample(S = 1, sim = sim, rr = r, tMax = 10)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, sim$TE[1]),
#'      xlab = "Mya")
#' lines(t, rev(r(t)))
#' 
#' ###
#' # now let us try a step function
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation was short lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we can create the sampling rate here from a few vectors
#' 
#' # rates
#' rlist <- c(1, 3, 0.5)
#' 
#' # rate shift times -  this could be c(10, 6, 2)
#' # and would produce the same function
#' rShifts <- c(0, 4, 8)
#' 
#' # create the rate to visualize it
#' r <- make.rate(rlist, tMax = 10, fShifts = rShifts)
#' 
#' # time
#' t <- seq(0, 10, by = 0.1)
#' 
#' # visualizing the plot from past to present
#' plot(x = t, y = rev(r(t)), main = "Simulated preservation", type = "l",
#'      xlab = "Mya", ylab = "preservation rate",
#'      xlim = c(10, sim$TE[1]))
#' 
#' # sample
#' occs <- sample(S = 1, sim = sim, rr = r, tMax = 10)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, sim$TE[1]),
#'      xlab = "Mya")
#' 
#' # frontiers of each regime
#' abline(v = 10 - rShifts, col = "red")
#' 
#' ###
#' # we can create a step function in a different way as well
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation was short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # preservation function
#' r <- function(t) {
#'   ifelse(t < 4, 1,
#'          ifelse(t < 8, 3, 0.5))
#' }
#' # note how this function should be exactly the same as the previous one
#' 
#' # time
#' t <- seq(0, 10, by = 0.1)
#' 
#' # visualizing the plot from past to present
#' plot(x = t, y = rev(r(t)), main = "Simulated preservation", type = "l",
#'      xlab = "Mya", ylab = "preservation rate",
#'      xlim = c(10, sim$TE[1]))
#' 
#' # sample
#' occs <- sample(S = 1, sim = sim, rr = r, tMax = 10)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, sim$TE[1]),
#'      xlab = "Mya")
#' abline(v = 10 - rShifts, col = "red")
#' 
#' ###
#' # finally we could generate sampling dependent on temperature
#' 
#' if (requireNamespace("RPANDA", quietly = TRUE)) {
#'   # simulate a group
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#'   
#'   # in case first simulation was short-lived
#'   while ((sim$TS[1] - sim$TE[1]) < 10) {
#'     sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#'   }
#'   
#'   # preservation function dependent on temperature
#'   r_t <- function(t, env) {
#'     return(0.25*env)
#'   }
#'   
#'   # get the temperature data
#'   data(InfTemp, package = "RPANDA")
#'   
#'   # final preservation
#'   r <- make.rate(r_t, envF = InfTemp)
#'   
#'   # visualizing the plot from past to present
#'   plot(x = t, y = rev(r(t)), main = "Simulated preservation", type = "l",
#'        xlab = "Mya", ylab = "preservation rate",
#'        xlim = c(10, sim$TE[1]))
#'   
#'   # sample
#'   occs <- sample(S = 1, sim = sim, rr = r, tMax = 10)
#'   
#'   # check histogram
#'   hist(occs,
#'        xlim = c(10, sim$TE[1]),
#'        xlab = "Mya")
#'   lines(t, rev(r(t)))
#' }
#' 
#' @name sample
#' @rdname sample
#' @export

sample<-function(S, sim, rr, tMax) {
  # invert times since simulation goes from 0 to tMax
  TE <- tMax - sim$TE
  TS <- tMax - sim$TS

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
    # take the waiting time for sampling, using rexp.var()
    WaitTimeR <- ifelse(rr(Now) > 0, rexp.var(1, rr, Now, tMax), Inf)

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
