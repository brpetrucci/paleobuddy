#' Non-constant rate Birth-Death simulation
#'
#' Simulates a species birth-death process with general rates for any number of
#' starting species. Allows for the speciation/extinction rate to be constant or
#' a function of time. Takes an optional shape argument for speciation and/or 
#' extinction, under which the corresponding rate will be taken to be a Weibull 
#' scale for age-dependent dynamics. Returns an object containing lists of 
#' speciation times, extinction times, parents and status (extant or not). Can 
#' return true extinction times or simply information on whether species lived
#' after maximum simulation time. Allows for constraining on the number of species
#' at the end of the simulation, either total or extant. Returns \code{NA} and 
#' sends a warning if it cannot find a simulation with the desired number of 
#' species after \code{100000} tries.
#' Note that while time runs from \code{0} to \code{tmax} on the function itself,
#' it runs from \code{tmax} to \code{0} on the lists returned to conform with the
#' literature. 
#'
#' @param n0 Initial number of species, usually 1. Good parameter
#' to tweak if one is observing a low sample size when testing.
#'
#' @param pp Function to hold the speciation rate over time. It will either be
#' interpreted as an exponential rate, or a Weibull scale if 
#' \code{pShape != NULL}.
#'
#' @param qq Similar to above, but for the extinction rate.
#' 
#' Note: this function is meant to be called by \code{bd.sim}, so it neither
#' allows for as much flexibility, nor calls \code{make.rate}. If the user wishes
#' to use \code{bd.sim.general} with environmental or step-function rates, they
#' can generate the rate with \code{make.rate} and supply it to the function.
#'
#' @param tMax Ending time of simulation. Any species still living after
#' \code{tMax} is considered extant, and any species that would be generated
#' after \code{tMax} is not born.
#'
#' @param pShape Shape parameter for the Weibull distribution for age-dependent
#' speciation. Default is \code{NULL}, where \code{pp} will be considered a
#' (possibly) time-dependent exponential rate. For \code{pShape != NULL}, 
#' \code{pp} will be considered a scale, and \code{rexp.var} will draw a Weibull 
#' distribution instead. 
#'
#' @param qShape similar to above, but for the extinction rate.
#' 
#' Note: Time-varying shape is implemented, so one could have \code{pShape} or
#' \code{qShape} be a function of time. It is not thoroughly tested, however, so 
#' it may be prudent to wait for a future release where this feature is well 
#' established.
#' 
#' @param nFinal An interval of acceptable number of species at the end of the
#' simulation. If not supplied, default is \code{c(0, Inf)}, so that any number
#' of species is accepted. If supplied, \code{bd.sim.general} will run until the
#' number of total species generated, or, if \code{extOnly = TRUE}, the number of
#' extant species at the end of the simulation, lies within the interval.
#' 
#' @param extOnly A boolean indicating whether \code{nFinal} should be taken as
#' the number of total or extant species during the simulation. If \code{TRUE},
#' \code{bd.sim.general} will run until the number of extant species lies within
#' the \code{nFinal} interval. If \code{FALSE}, as default, it will run until the
#' total number of species generated lies within that interval.
#' 
#' @param fast When \code{TRUE}, sets \code{rexp.var} to throw away waiting times 
#' higher than the maximum simulation time. Should be \code{FALSE} for unbiased 
#' testing of age dependency. User might also se it to \code{FALSE} for more 
#' accurate waiting times.
#' 
#' @param trueExt When \code{TRUE}, time of extinction of extant species will be 
#' the true time, otherwise it will be \code{tMax+0.01}. Need to be \code{TRUE} 
#' when testing age-dependent extinction.
#'
#' @return A list of vectors, as follows
#'
#' \describe{
#' \item{\code{TE}}{List of extinction times, with \code{-0.01} as the time of
#' extinction for extant species.}
#'
#' \item{\code{TS}}{List of speciation times, with \code{tMax + 0.01} as the time 
#' of speciation for species that started the simulation.}
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
#' # we can test a couple scenarios
#' 
#' ###
#' # first, even though this is bd.sim.general, we can try constant rates
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- 0.11
#' 
#' # extinction
#' q <- 0.08
#' 
#' # run the simulation
#' sim <- bd.sim.general(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # we can complicate things further with a linear function as a rate
#' # bd.sim.general takes longer so we run examples for 1000 replicates instead
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- function(t) {
#'   return(0.03 + 0.005*t)
#' }
#' 
#' # extinction
#' q <- 0.05
#' 
#' # run the simulation
#' sim <- bd.sim.general(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # we can also create a step function
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation rate
#' p <- function(t) {
#'   return(0.03 + 0.005*t)
#' }
#' 
#' # list of extinction rates
#' qList <- c(0.06, 0.09, 0.11)
#' 
#' # list of shift times. Note qShifts could be c(40, 20, 10) for
#' # identical results
#' qShifts <- c(0, 15, 25)
#' 
#' # let us take a look at how make.rate will make it a step function
#' q <- make.rate(qList, tMax = tMax, fShifts = qShifts)
#' 
#' # and plot it
#' plot(seq(0, tMax, 0.1), q(seq(0, tMax, 0.1)), type = 'l',
#'      main = "Extintion rate as a step function", xlab = "Time (My)",
#'      ylab = "Rate (species/My)")
#' 
#' # a different way to define the same extinction function
#' q <- function(t) {
#'   ifelse(t < 15, 0.06,
#'          ifelse(t < 25, 0.09, 0.11))
#' }
#' 
#' # run the simulation
#' sim <- bd.sim.general(n0, p, q, tMax, nFinal = c(2, Inf))
#' # we could instead have used q made with make.rate
#' # that is, however, much slower
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # another feature to add is age dependency
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- 0.1
#' 
#' # extinction - a Weibull scale
#' q <- 10
#' 
#' # extinction shape
#' qShape <- 1
#' 
#' # run simulations
#' sim <- bd.sim.general(n0, p, q, tMax, qShape = qShape, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # scale can be time-dependent
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- 0.15
#' 
#' # extinction - a Weibull scale
#' q <- function(t) {
#'   return(8 + 0.05*t)
#' }
#' 
#' # extinction shape
#' qShape <- 1
#' 
#' # run simulations
#' sim <- bd.sim.general(n0, p, q, tMax, qShape = qShape, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # finally, we could have environmental dependency on a rate
#' if (requireNamespace("RPANDA", quietly = TRUE)) {
#'   # initial number of species
#'  n0 <- 1
#'
#'  # maximum simulation time
#'  tMax <- 40
#'  
#'  # temperature-dependent speciation
#'  p_t <- function(t, temp) {
#'    return(0.04*exp(0.15*temp))
#'  }
#'
#'  # extinction
#'  q <- 0.075
#'
#'  # using RPANDA to get the temperature data
#'  data(InfTemp, package="RPANDA")
#'
#'  # speciation
#'  p <- make.rate(p_t, envF = InfTemp)
#'
#'  # run simulations
#'  sim <- bd.sim.general(n0, p, q, tMax, nFinal = c(2, Inf))
#'   
#'   # we can plot the phylogeny to take a look
#'   if (requireNamespace("ape", quietly = TRUE)) {
#'     phy <- make.phylo(sim)
#'     ape::plot.phylo(phy)
#'   }
#' }
#' 
#' # note nFinal has to be sensible
#' \dontrun{
#' # this would return a warning, since it is virtually impossible to get 100
#' # species at a process with diversification rate -0.09 starting at n0 = 1
#' sim <- bd.sim.general(1, pp = 0.01, qq = 1, tMax = 100, nFinal = c(100, Inf))
#' }
#'
#' @name bd.sim.general
#' @rdname bd.sim.general
#' @export
#' 

bd.sim.general <- function(n0, pp, qq, tMax, 
                         pShape = NULL, qShape = NULL,
                         nFinal = c(0, Inf), extOnly = FALSE,
                         fast = TRUE, trueExt = FALSE) {
  # initialize species count with a value that makes sure the while loop runs
  len <- -1
  
  # counter to make sure the nFinal is achievable
  counter <- 1
  
  # if shape is not null, make scale a function to facilitate checking
  if (!is.null(pShape) & is.numeric(pp)) {
    p <- pp
    pp <- Vectorize(function(t) p)
  }
  
  if (!is.null(qShape) & is.numeric(qq)) {
    q <- qq
    qq <- Vectorize(function(t) q)
  }
  
  while (len < nFinal[1] | len > nFinal[2]) {
    # create vectors to hold times of speciation, extinction, parents and status
    TS <- rep(-0.01, n0)
    TE <- rep(NA, n0)
    parent <- rep(NA, n0)
    isExtant <- rep(TRUE, n0)
  
    # initialize species count
    sCount <- 1
  
    # while we have species to be analyzed still
    while (length(TE) >= sCount) {
      # actual speciation time
      specT <- ifelse(TS[sCount] < 0, 0, TS[sCount])
      
      # need this to pass speciation time to rexp.var only when shape is there
      pSpecT <- NULL
      if (!is.null(pShape)) pSpecT <- specT
      
      qSpecT <- NULL
      if (!is.null(qShape)) qSpecT <- specT
  
      # get the time of speciation, or 0 if the species
      # was there at the beginning
      tNow <- specT

      # find the waiting time using rexp.var if pp is not constant
      # note we need to pass NULL for TS if the corresponding shape is NULL
      waitTimeS <- ifelse(
        is.numeric(pp), rexp(1, pp), 
        ifelse(pp(tNow) > 0, 
               rexp.var(1, pp, tNow, tMax, pShape, pSpecT, fast), Inf))
      waitTimeE <- ifelse(
        is.numeric(qq), rexp(1, qq), 
        ifelse(qq(tNow) > 0,
               rexp.var(1, qq, tNow, tMax, qShape, qSpecT, fast), Inf))
  
      tExp <- tNow + waitTimeE
  
      # while there are fast enough speciations before the species goes extinct,
      while ((tNow + waitTimeS) <= min(tExp, tMax)) {
  
        # advance to the time of speciation
        tNow <- tNow + waitTimeS
  
        # add new times to the vectors
        TS <- c(TS, tNow)
        TE <- c(TE, NA)
        parent <- c(parent, sCount)
        isExtant <- c(isExtant, TRUE)
  
        # get a new speciation waiting time, and include it in the vector
        waitTimeS <- ifelse(
          is.numeric(pp), rexp(1, pp), 
          ifelse(pp(tNow) > 0, 
                 rexp.var(1, pp, tNow, tMax, pShape, pSpecT, fast), Inf))
      }
  
      # reached the time of extinction
      tNow <- tExp
  
      # record extinction, and if species is extant make it more than tMax
      TE[sCount] <- ifelse(tNow < tMax | trueExt, tNow, tMax + 0.01)
      isExtant[sCount] <- ifelse(TE[sCount] > tMax, TRUE, FALSE)
  
      # next species
      sCount <- sCount + 1
    }
  
    # now we invert TE and TS so time goes from tMax to 0
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
