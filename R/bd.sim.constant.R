#' Constant rate Birth-Death simulation
#'
#' Simulates a species birth-death process with constant rates for any number
#' of starting species. Returns an object containing lists of speciation times,
#' extinction times, parents and status (extant or not). Allows for constraining
#' on the number of species at the end of the simulation, either total or extant.
#' Note that while time runs from \code{0} to \code{tmax} on the function itself,
#' it runs from \code{tmax} to \code{0} on the lists returned to conform with the
#' literature. Returns \code{NA} and sends a warning if it cannot find a
#' simulation with the desired number of species after \code{100000} tries.
#'
#' @param n0 Initial number of species, usually 1. Good parameter
#' to tweak if one is observing a low sample size when testing.
#'
#' @param pp Speciation rate. Must be constant.
#'
#' @param qq Extinction rate, similar to above.
#' 
#' Note: \code{pp} and \code{qq} must always be greater than 0
#'
#' @param tMax Ending time of simulation. Any species still living after
#' \code{tMax} is considered extant, and any species that would be generated
#' after \code{tMax} is not born.
#' 
#' @param nFinal An interval of acceptable number of species at the end of the
#' simulation. If not supplied, default is \code{c(0, Inf)}, so that any number
#' of species is accepted. If supplied, \code{bd.sim.constant} will run until the
#' number of total species generated, or, if \code{extOnly = TRUE}, the number of
#' extant species at the end of the simulation, lies within the interval.
#' 
#' @param extOnly A boolean indicating whether \code{nFinal} should be taken as
#' the number of total or extant species during the simulation. If \code{TRUE},
#' \code{bd.sim.constant} will run until the number of extant species lies within
#' the \code{nFinal} interval. If \code{FALSE}, as default, it will run until the
#' total number of species generated lies within that interval.
#'
#' @return A list of vectors, as follows
#'
#' \describe{
#' \item{\code{TE}}{List of extinction times, with -0.01 as the time of
#' extinction for extant species.}
#'
#' \item{\code{TS}}{List of speciation times, with tMax+0.01 as the time of
#' speciation for species that started the simulation.}
#'
#' \item{\code{PAR}}{List of parents. Species that started the simulation have
#' NA, while species that were generated during the simulation have their
#' parent's number. Species are numbered as they are born.}
#'
#' \item{\code{EXTANT}}{List of booleans representing whether each species is
#' extant.}}
#'
#' @author Bruno do Rosario Petrucci.
#'
#' @examples
#' 
#' # we can show a couple scenarios
#' 
#' ###
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
#' # run the simulation
#' sim <- bd.sim.constant(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
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
#' # run the simulation
#' sim <- bd.sim.constant(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
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
#' # run the simulation
#' sim <- bd.sim.constant(n0, p, q, tMax)
#' 
#' # of course in this case there are no phylogenies to plot
#' 
#' # note nFinal has to be sensible
#' \dontrun{
#' # this would return an error
#' sim <- bd.sim.constant(1, pp = 0.01, qq = 1, tMax = 100, nFinal = c(100, Inf))
#' }
#' 
#' @name bd.sim.constant
#' @rdname bd.sim.constant
#' @export


bd.sim.constant <- function(n0 = 1, pp, qq, tMax, 
                            nFinal = c(0, Inf), extOnly = FALSE) {
  # initialize species count with a value that makes sure the while loop runs
  len <- -1
  
  # counter to make sure the nFinal is achievable
  counter <- 1
  
  while (len < nFinal[1] | len > nFinal[2]) {
    # check that the rates are constant
    if (!(is.numeric(pp) & length(pp) == 1 &
        is.numeric(qq) * length(qq) == 1)) {
      stop("bd.sim.constant requires constant rates")
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
  
    # check the size of the simulation
    len <- ifelse(extOnly, sum(isExtant), length(isExtant))
    # if this is in nFinal, the while loop stops
    
    # if we have ran for too long, stop
    counter <- counter + 1
    if (counter > 100000) {
      warning("This value of nFinal took more than 100000 simulations 
              to achieve")
      return(NA)
    }
  }

  return(list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant))
}
