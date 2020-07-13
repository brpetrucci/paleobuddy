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
#' @param nFinal an interval of acceptable number of species at the end of the
#' simulation. If not supplied, default is \code{c(0, Inf)}, so that any number
#' of species is accepted. If supplied, \code{BDSimConstant} will run until the
#' number of total species generated, or, if \code{extOnly = TRUE}, the number of
#' extant species at the end of the simulation, lies within the interval.
#' 
#' @param extOnly a boolean indicating whether \code{nFinal} should be taken as
#' the number of total or extant species during the simulation. If \code{TRUE},
#' \code{BDSimConstant} will run until the number of extant species lies within
#' the \code{nFinal} interval. If \code{FALSE}, as default, it will run until the
#' total number of species generated lies within that interval.
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
#' # first we define a time vector for plotting
#' tMax <- 40
#' time <- 1:tMax
#' 
#' # now we can show a couple scenarios
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
#' sim <- BDSim(n0, p, q, tMax)
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
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
#' sim <- BDSimConstant(n0, p, q, tMax)
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
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
#' sim <- BDSim(n0, p, q, tMax)
#' 
#' # of course in this case there are no phylogenies to plot
#' 

#'
#' @name BDSimConstant
#' @rdname BDSimConstant
#' @export


BDSimConstant <- function(n0 = 1, pp, qq, tMax, 
                          nFinal = c(0, Inf), extOnly = FALSE) {
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
  
  # check if we do not have an acceptable number of species
  test <- (extOnly & (sum(isExtant) < nFinal[1] | sum(isExtant) > nFinal[2])) |
    (length(TE) < nFinal[1] | length(TE) > nFinal[2])
  
  # if we do not, find a new result
  if (test) {
    return(BDSimConstant(n0, pp, qq, tMax, nFinal, extOnly))
  }

  # if we do, return this one
  else {
    return(list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant))
  }
}