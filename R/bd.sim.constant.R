#' Constant rate Birth-Death simulation
#' 
#' Simulates a species birth-death process with general rates for any number of
#' starting species. Allows for constraining results on the number of species at 
#' the end of the simulation, either total or extant. Returns an object containing 
#' vectors of speciation times, extinction times, parents (= species' mother 
#' species) and status at the end of the simulation (extant or not) for each 
#' species in the simulation. It may return true extinction times or simply 
#' information on whether species lived after the maximum simulation time. For 
#' time-varying and age-dependent simulations, see \code{bd.sim.general}. For a 
#' function that unites both scenarios, see \code{bd.sim}.
#' Please note while time runs from \code{0} to \code{tMax} in the simulation, it 
#' returns speciation/extinction times as \code{tMax} (origin of the group) to 
#' \code{0} (the "present" and end of simulation), so as to conform to other
#' packages in the literature.
#' 
#' @inheritParams bd.sim 
#' 
#' @param pp Speciation rate. Must be constant.
#'
#' @param qq Extinction rate, similar to above.
#' 
#' Note: \code{pp} and \code{qq} must always be greater than 0
#'
#' @return A list of vectors, as follows
#'
#' \describe{
#' \item{\code{TE}}{List of extinction times, with \code{-0.01} as the time of
#' extinction for extant species.}
#'
#' \item{\code{TS}}{List of speciation times, with \code{tMax+0.01} as the time of
#' speciation for species that started the simulation.}
#'
#' \item{\code{PAR}}{List of parents. Species that started the simulation have
#' \code{NA}, while species that were generated during the simulation have their
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
#' # this would return a warning, since it is virtually impossible to get 100
#' # species at a process with diversification rate -0.09 starting at n0 = 1
#' sim <- bd.sim.constant(1, pp = 0.01, qq = 1, tMax = 100, nFinal = c(100, Inf))
#' }
#' 
#' @name bd.sim.constant
#' @rdname bd.sim.constant
#' @export


bd.sim.constant <- function(n0, pp, qq, tMax, 
                            nFinal = c(0, Inf), extOnly = FALSE,
                            trueExt = FALSE) {
  # check that the rates are constant
  if (!(is.numeric(pp) & length(pp) == 1 &
        is.numeric(qq) * length(qq) == 1)) {
    stop("bd.sim.constant requires constant rates")
  }
  
  # check that the rates are non-negative
  if (pp < 0 || qq < 0) {
    stop("rates cannot be negative")
  }
  
  # check that n0 is not negative
  if (n0 <= 0) {
    stop("initial number of species must be positive")
  }
  
  # check nFinal's length
  if (length(nFinal) != 2) {
    stop("nFinal must be a vector with a minimum and maximum number 
         of species")
  }
  
  # initialize species count with a value that makes sure the while loop runs
  len <- -1
  
  # counter to make sure the nFinal is achievable
  counter <- 1
  
  while (len < nFinal[1] | len > nFinal[2]) {
    # if we have gone through more than 100000 simulations, it is probably 
    # impossible to achieve the nFinal we want
    if (counter > 100000) {
      warning("This value of nFinal took more than 100000 simulations 
              to achieve")
      return(NA)
    }
    
    # initialize the vectors to hold times of speciation and extinction, parents
    # and status (extant or not)
    TS <- rep(0, n0)
    TE <- rep(NA, n0)
    parent <- rep(NA, n0)
    isExtant <- rep(TRUE, n0)
  
    # initialize the counting variable
    sCount <- 1
  
    # while we have more species in a vector than we have analyzed,
    while (length(TE) >= sCount) {
      # start at the time of speciation of sCount
      tNow <- TS[sCount]
  
      # draw waiting times with rexp()
      waitTimeS <- ifelse(pp > 0, rexp(1, pp), Inf)
      waitTimeE <- ifelse(qq > 0, rexp(1, qq), Inf)
  
      # if the time of extinction is after the end of the simulation, make it tMax
      tExp <- tNow + waitTimeE
  
      # while there are fast enough speciations before the species goes extinct,
      while ((tNow + waitTimeS) <= min(tExp, tMax)) {
        # update time
        tNow <- tNow + waitTimeS
  
        # create a new species with corresponding TE, TS and parent
        TS <- c(TS, tNow)
        TE <- c(TE, NA)
        parent <- c(parent, sCount)
        isExtant <- c(isExtant, TRUE) # it is alive
  
        # take a new waiting time - if now + waitTimeS is still less than when
        # the species goes extinct, repeat
        waitTimeS <- ifelse(pp > 0, rexp(1, pp), Inf)
      }
  
      # reached the time of the species extinction
      tNow <- tExp
  
      # record the extinction - if TE[sCount] is NA, it is extant
      isExtant[sCount] <- ifelse(is.na(TE[sCount]), TRUE, FALSE)
  
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
  }

  return(list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant))
}
