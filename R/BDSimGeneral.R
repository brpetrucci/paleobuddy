#' Non-constant rate Birth-Death simulation
#'
#' \code{BDSimGeneral} takes an initial number of species, speciation and
#' extinction rates (either constants, functions of time, or of time and an
#' environmental variable), a maximum simulation time and possibly a shape for
#' age-dependent speciation and/or extinction. It then generates the speciation 
#' and extinction times, parent and status information for the species in the 
#' simulation. Time runs from \code{tMax} to 0, to be consistent with the 
#' literature, though one can easily invert that by subtracting the results from
#' \code{tMax}.
#'
#' @param n0 initial number of species, usually 1. Good parameter
#' to tweak if one is observing a low sample size when testing.
#'
#' @param pp function to hold the speciation rate over time. It will either be
#' interpreted as an exponential rate, or a Weibull scale if 
#' \code{pShape != NULL}.
#'
#' @param qq similar to above, but for the extinction rate.
#' 
#' Note: this function is meant to be called by \code{BDSim}, so it neither
#' allows for as much flexibility, nor call \code{MakeRate}. If the user wishes
#' to use \code{BDSimGeneral} with environmental or step-function rates, they
#' can generate the rate with \code{MakeRate} and supply it to the function.
#'
#' @param tMax ending time of simulation. Any species still living
#' after \code{tMax} is considered extant, and any species that would be
#' generated after \code{tMax} is not born.
#'
#' @param pShape shape parameter for the Weibull distribution for age-dependent
#' speciation. Default is \code{NULL}, where \code{pp} will be considered a
#' time-dependent exponential rate. For \code{pShape != NULL}, \code{pp} will
#' be considered a scale, and \code{rexp_var} will draw a Weibull distribution
#' instead.
#'
#' @param qShape similar as above, but for the extinction rate.
#' 
#' @param nFinal an interval of acceptable number of species at the end of the
#' simulation. If not supplied, default is \code{c(0, Inf)}, so that any number
#' of species is accepted. If supplied, \code{BDSimGeneral} will run until the
#' number of total species generated, or, if \code{extOnly = TRUE}, the number of
#' extant species at the end of the simulation, lies within the interval.
#' 
#' @param extOnly a boolean indicating whether \code{nFinal} should be taken as
#' the number of total or extant species during the simulation. If \code{TRUE},
#' \code{BDSimGeneral} will run until the number of extant species lies within
#' the \code{nFinal} interval. If \code{FALSE}, as default, it will run until the
#' total number of species generated lies within that interval.
#' 
#' @param fast used for \code{BDSimGeneral}. When \code{TRUE}, sets 
#' \code{rexp_var} to throw away waiting times higher than the maximum 
#' simulation time. Should be \code{FALSE} for unbiased testing of age 
#' dependency. User might also se it to \code{FALSE} for more accurate waiting
#' times.
#' 
#' @param trueExt used for \code{BDSimGeneral}. When \code{TRUE}, time of 
#' extinction of extant species will be the true time, otherwise it will be 
#' tMax+0.01. Need to be \code{TRUE} when testing age-dependent 
#' extinction.
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
#' # we can test a couple scenarios
#' 
#' ###
#' # first, even though this is BDSimGeneral, we can try constant rates
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
#' sim <- BDSimGeneral(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # we can complicate things further with a linear function as a rate
#' # BDSimGeneral takes longer so we run examples for 1000 replicates instead
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- function(t) {
#'   return(0.05 + 0.005*t)
#' }
#' 
#' # extinction
#' q <- 0.05
#' 
#' # run the simulation
#' sim <- BDSimGeneral(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # we can also create a step function. Keep in mind this is a slower way than by
#' # creating step functions using ifelse()
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
#' qList <- c(0.05, 0.08, 0.11)
#' 
#' # list of shift times. Note qShifts could be c(40, 20, 10) for
#' # identical results
#' qShifts <- c(0, 15, 25)
#' 
#' # let us take a look at how MakeRate will make it a step function
#' q <- MakeRate(qList, fShifts = qShifts)
#' 
#' # and plot it
#' plot(seq(0, tMax, 0.1), q(seq(0, tMax, 0.1)), type = 'l',
#'      main = "Extintion rate as a step function", xlab = "Time (My)",
#'      ylab = "Rate (species/My)")
#' 
#' # a different way to define the same extinction function
#' q <- function(t) {
#'   ifelse(t < 15, 0.05,
#'          ifelse(t < 25, 0.08, 0.11))
#' }
#' 
#' # run the simulation
#' sim <- BDSimGeneral(n0, p, q, tMax, nFinal = c(2, Inf))
#' # equivalent:
#' # sim <- BDSimGeneral(n0, p, qList, tMax, qShifts = qShifts)
#' # this is, however, much slower
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
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
#' sim <- BDSimGeneral(n0, p, q, tMax, qShape = qShape, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
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
#' sim <- BDSimGeneral(n0, p, q, tMax, qShape = qShape, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # finally, we could have environmental dependency on a rate
#' if (requireNamespace("RPANDA", quietly = TRUE)) {
#'   # initial number of species
#'   n0 <- 1
#'   
#'   # maximum simulation time
#'   tMax <- 40
#'   
#'   # temperature-dependent speciation
#'   p_t <- function(t, temp) {
#'     return(0.04*exp(0.15*temp))
#'   }
#'   
#'   # extinction
#'   q <- 0.075
#'   
#'   # using RPANDA to get the temperature data
#'   data(InfTemp, package="RPANDA")
#'   
#'   # speciation
#'   p <- MakeRate(p_t, tMax, envF = InfTemp)
#'   
#'   # since we need many species to be able to test this effectively using
#'   # RPANDA, and the rates become really noisy with temperature, we set
#'   # only 100 simulations to finish it in a reasonable time
#'   
#'   # run simulations
#'   sim <- BDSimGeneral(n0, p, q, tMax, nFinal = c(2, Inf))
#'   
#'   # we can plot the phylogeny to take a look
#'   if (requireNamespace("ape", quietly = TRUE)) {
#'     phy <- MakePhylo(sim)
#'     ape::plot.phylo(phy)
#'   }
#' }
#'
#' @name BDSimGeneral
#' @rdname BDSimGeneral
#' @export
#' 

BDSimGeneral <- function(n0, pp, qq, tMax, 
                         pShape = NULL, qShape = NULL,
                         nFinal = c(0, Inf), extOnly = FALSE,
                         fast = TRUE, trueExt = FALSE) {
  # initialize species count with a value that makes sure the while loop runs
  len <- -1
  
  # counter to make sure the nFinal is achievable
  counter <- 1
  
  # if shape is not null, make scale a function to facilitate checking
  if (!is.null(pShape)) {
    p <- pp
    pp <- ifelse(is.numeric(p), Vectorize(function(t) p), p)
  }
  
  if (!is.null(qShape)) {
    q <- qq
    qq <- ifelse(is.numeric(q), Vectorize(function(t) q), q)
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
  
      # get the time of speciation, or 0 if the species
      # was there at the beginning
      tNow <- ifelse(TS[sCount] < 0, 0, TS[sCount])

      # find the waiting time using rexp_var - note that in rexp_var we only
      # count t from tNow (to consider the rates as functions), so that
      # now we need to subtract tNow
      waitTimeS <- ifelse(
        is.numeric(pp), rexp(1, pp), 
        ifelse(pp(tNow) > 0, 
               rexp_var(1, pp, tNow, tMax, pShape, 
                        ifelse(TS[sCount] < 0, 0, TS[sCount]), fast), Inf))
      
      waitTimeE <- ifelse(
        is.numeric(qq), rexp(1, qq), 
        ifelse(qq(tNow) > 0,
               rexp_var(1, qq, tNow, tMax, qShape,
                        ifelse(TS[sCount] < 0, 0, TS[sCount]), fast), Inf))
  
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
                 rexp_var(1, pp, tNow, tMax, pShape, 
                          ifelse(TS[sCount] < 0, 0, TS[sCount]), fast), Inf))
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
      message("This value of nFinal took more than 100000 simulations 
              to achieve")
      return(NA)
    }
  }
  
  return(list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant))
}
