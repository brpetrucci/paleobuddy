#' Constant rate Birth-Death simulation
#'
#' \code{BDSimConstant} takes a rate, a maximum simulation time and an initial
#' number of species and returns the history of the clade originated from those
#' species using a birth-death process with constant rates equal to the given 
#' rates. It then generates the speciation and extinction times, parent and 
#' status information for the species in the simulation. Time runs from 
#' \code{tMax} to 0, to be consistent with the literature, though one can easily 
#' invert that by subtracting the results from \code{tMax}.
#'
#' @param n0 initial number of species, usually 1. Good parameter
#' to tweak if one is observing a low sample size when testing.
#'
#' @param pp speciation rate. Must be constant.
#'
#' @param qq extinction rate, similar to above.
#' 
#' Note: \code{pp} and \code{qq} must always be greater than 0
#'
#' @param tMax ending time of simulation. Any species still living
#' after tMax is considered extant, and any species that would be generated
#' after \code{tMax} is not born.
#'
#' @return a list of vectors, as follows
#'
#' \describe{
#' \item{\code{TE}}{list of extinction times, with -0.01 as the time of
#' extinction for extant species.}
#'
#' \item{\code{TS}}{list of speciation times, with tMax+0.01 as the time of
#' speciation for species that started the simulation.}
#'
#' \item{\code{PAR}}{list of parents. Species that started the simulation have
#' NA, while species that were generated during the simulation have their
#' parent's number. Species are numbered as they are born.}
#'
#' \item{\code{EXTANT}}{list of booleans representing whether each species is
#' extant.}}
#'
#' @author written by Bruno do Rosario Petrucci.
#'
#' @examples
#' 
#' # first we define a function to calculate the mean diversity and var at time t
#' sim.mean<-function(t, simList) {
#'   # apply to each replicate of the simulation
#'   simExtantT<-as.numeric(lapply(1:length(simList), function(y) {
#'     
#'     # invert extinction and speciation times, so that we may consider functions
#'     # running from past to present
#'     TS <- tMax - simList[[y]]$TS
#'     TE <- tMax - simList[[y]]$TE
#'     
#'     # find how many species are alive at time t
#'     length(which(TS <= t & TE >= t))}))
#'   return(list(mean = mean(simExtantT), var = (var(simExtantT))))
#' }
#' 
#' # also, we need functions to calculate the expected variance at time t
#' int<-function(t, div) {
#'   # integrate diversity from 0 to time t
#'   return(integrate(div, 0, t)$value)
#' }
#' # formula for the variance - see Kendall 1948
#' div.var<-function(t, div, qq) {
#'   return(n0*exp(int(t, div))*(exp(int(t, div)) - 1 + 2*exp(int(t, div))*
#'                                 integrate(Vectorize(function(x)
#'                                   exp(-int(x, div))*qq(x)), 0, t)$value))
#' }
#' 
#' # and a time parameter we will need for plotting
#' tMax <- 40
#' time <- 1:tMax
#' 
#' # now we can test a couple scenarios
#' # first, extinction 0
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- 0.1
#' 
#' # extinction
#' q <- 0
#' 
#' \dontrun{ 
#' # run the simulations
#' simList <- lapply(1:10000, function(x) BDSimConstant(n0, p, q, tMax))
#' 
#' 
#' # let us make vectors to hold the average diversity and variance
#' 
#' # function for diversity
#' div <- Vectorize(function(t) p - q)
#' 
#' # function for extinction - needed for the variance
#' qq <- Vectorize(function(t) q)
#' 
#' # calculate the mean diversity at our time points
#' meanDiv <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$mean))
#' 
#' # calculate the expected diversity for the sime time points
#' expectedDiv <- VarRateExp(div, 1, time)
#' 
#' # do the same with variance
#' meanVar <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$var))
#' expectedVar <- unlist(lapply(time, function(x) div.var(x, div, qq)))
#' 
#' # and now let us check out the plots
#' 
#' # mean diversity, compared with expected
#' plot(time, log(meanDiv), type = 'l', main = "Species diversity", 
#'      xlab = "Time (My)", ylab = "log(Diversity)")
#' lines(time, log(expectedDiv), col = 'RED')
#' legend(x = 5, y = log(max(meanDiv)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' 
#' # same for variance
#' plot(time, log(meanVar), type = 'l',  main = "Species diversity variance", 
#'      xlab = "Time (My)", ylab = "log(Variance)")
#' lines(time, log(expectedVar), type = 'l', col='RED')
#' legend(x = 5, y = log(max(meanVar)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' }
#' 
#' # now let us try to turn extinction up a bit
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- 0.1
#' 
#' # extinction
#' q <- 0.04
#' 
#' \dontrun{
#' # run the simulations
#' simList <- lapply(1:10000, function(x) BDSimConstant(n0, p, q, tMax))
#' 
#' # calculate the mean diversity at our time points
#' meanDiv <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$mean))
#' 
#' # calculate the expected diversity for the sime time points
#' expectedDiv <- VarRateExp(div, 1, time)
#' 
#' # do the same with variance
#' meanVar <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$var))
#' expectedVar <- unlist(lapply(time, function(x) div.var(x, div, qq)))
#' 
#' # and now let us check out the plots
#' 
#' # mean diversity, compared with expected
#' plot(time, log(meanDiv), type = 'l', main = "Species diversity", 
#'      xlab = "Time (My)", ylab = "log(Diversity)")
#' lines(time, log(expectedDiv), col = 'RED')
#' legend(x = 5, y = log(max(meanDiv)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' 
#' # same for variance
#' plot(time, log(meanVar), type = 'l',  main = "Species diversity variance", 
#'      xlab = "Time (My)", ylab = "log(Variance)")
#' lines(time, log(expectedVar), type = 'l', col='RED')
#' legend(x = 5, y = log(max(meanVar)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' }
#' 
#' # we can also try a pure-death process
#' 
#' # initial number of species - note the high number, so we get an appreciable
#' # sample size
#' n0 <- 100
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- 0
#' 
#' # extinction
#' q <- 0.02
#' 
#' \dontrun{
#' # run the simulations
#' simList <- lapply(1:10000, function(x) BDSimConstant(n0, p, q, tMax))
#' 
#' # calculate the mean diversity at our time points
#' meanDiv <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$mean))
#' 
#' # calculate the expected diversity for the sime time points
#' expectedDiv <- VarRateExp(div, 1, time)
#' 
#' # do the same with variance
#' meanVar <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$var))
#' expectedVar <- unlist(lapply(time, function(x) div.var(x, div, qq)))
#' 
#' # and now let us check out the plots
#' 
#' # mean diversity, compared with expected
#' plot(time, log(meanDiv), type = 'l', main = "Species diversity", 
#'      xlab = "Time (My)", ylab = "log(Diversity)")
#' lines(time, log(expectedDiv), col = 'RED')
#' legend(x = 30, y = log(max(meanDiv)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' 
#' # same for variance
#' plot(time, log(meanVar), type = 'l',  main = "Species diversity variance", 
#'      xlab = "Time (My)", ylab = "log(Variance)")
#' lines(time, log(expectedVar), type = 'l', col='RED')
#' legend(x = 30, y = log(max(meanVar)) - 2, legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' }
#' # all the cases seem to agree pretty well with expectations
#'
#' @name BDSimConstant
#' @rdname BDSimConstant
#' @export


BDSimConstant <- function(n0 = 1, pp, qq, tMax) {
  # check that the rates are constant
  if (!(is.numeric(pp) & length(pp) == 1 &
      is.numeric(qq) * length(qq) == 1)) {
    stop("BDSimConstant requires constant rates")
  }
  
  # initialize the vectors to hold times of speciation and extinction, parents
  # and status (extant or not)
  TS <- rep(-0.01, n0)
  TE <- rep(NA, n0)
  parent <- rep(NA, n0)
  isExtant <- rep(TRUE, n0)

  # initialize the counting variable
  sCount <- 1

  # while we have more species in a vector than we have analyzed,
  while (length(TE) >= sCount) {
    # TS starts at -0.01 to show it was alive at the beginning, but to count
    # time we need to start at 0
    tNow <- ifelse(TS[sCount] < 0, 0, TS[sCount])

    # draw waiting times with rexp()
    waittimeS <- ifelse(pp > 0, rexp(1, pp), Inf)
    waittimeE <- ifelse(qq > 0, rexp(1, qq), Inf)

    # if the time of extinction is after the end of the simulation, make it tMax
    tExp <- min(tNow + waittimeE, tMax)

    # while there are fast enough speciations before the species goes extinct,
    while ((tNow + waittimeS) <= tExp) {
      # update time
      tNow<-tNow + waittimeS

      # create a new species with corresponding TE, TS and parent
      TS <- c(TS, tNow)
      TE <- c(TE, NA)
      parent <- c(parent, sCount)
      isExtant <- c(isExtant, TRUE) # it is alive

      # take a new waiting time - if now + waittimeS is still less than when
      # the species goes extinct, repeat
      waittimeS <- ifelse(pp > 0, rexp(1, pp), Inf)
    }

    # reached the time of the species extinction
    tNow <- tExp

    # record the extinction - if tExp >= tMax, it didn't go extinct
    TE[sCount] <- ifelse(tNow < tMax, tNow, tMax + 0.01)
    isExtant[sCount] <- ifelse(TE[sCount] > tMax, TRUE, FALSE)

    # next species
    sCount <- sCount + 1
  }

  # finally, we invert both TE and TS to attain to the convention that time
  # runs from 0 to tMax
  TE <- tMax - TE
  TS <- tMax - TS

  return(list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant))
}

