#' General rate Birth-Death simulation
#'
#' \code{BDSim} takes an initial number of species, speciation and extinction
#' rate functions and a maximum time of simulation, together with multiple
#' options to alter the rates, and calls \code{BDSimConstant} or
#' \code{BDSimGeneral} to generate a species diversification process under the
#' desired scenario.
#'
#' @param n0 initial number of species, usually 1. Good parameter
#' to tweak if one is observing a low sample size when testing.
#'
#' @param pp function to hold the speciation rate over time. It could be a
#' constant or a function of time (to be an exponential rate or weibull scale),
#' a function of time and an environmental variable, or a vector of rates to be
#' accompanied by a vector of rate shifts \code{pShifts}.
#'
#' @param qq similar as above, but for the extinction rate.
#'
#' Note: \code{pp} and \code{qq} must always be greater than 0
#'
#' @param tMax ending time of simulation. Any species still living
#' after tMax is considered extant, and any species that would be generated
#' after \code{tMax} is not born.
#'
#' @param pShape shape parameter for the Weibull distribution for age-dependent
#' speciation. Default is \code{NULL}, where \code{pp} will be considered a
#' time-dependent exponential rate. For \code{pShape != NULL}, \code{pp} will
#' be considered a scale, and \code{rexp_var} will draw a Weibull distribution
#' instead.
#'
#' @param qShape similar as above, but for the extinction rate.
#'
#' @param envPP a matrix containing time points and values of an
#' environmental variable, like temperature, for each time point. This will be
#' used to create a speciation rate, so \code{pp} must be a function of time
#' and said variable.
#'
#' @param envQQ similar as above, but for the extinction rate.
#'
#' @param pShifts vector of rate shifts. First element must be the sstarting
#' time for the simulation (0 or tMax). It must have the same length as
#' \code{pp}. E.g. \code{pp = c(0.1, 0.2, 0.1)}, \code{pShifts = c(0, 10, 20)}
#' means the speciation rate will be 0.1 from 0 to 10, 0.2 from 10 to 20, and 
#' 0.1 from 20 to \code{tMax}. It would also be identical, in this case, to use
#' \code{pShifts = c(tMax, tMax-10, tMax-20)}.
#' 
#' Note that using this  method for step-function rates is currently slower than using
#' \code{ifelse}.
#'
#' @param qShifts similar as above, but for the extinction rate.
#' 
#' @param nFinal an interval of acceptable number of species at the end of the
#' simulation. If not supplied, default is \code{c(0, Inf)}, so that any number
#' of species is accepted. If supplied, \code{BDSimConstant} or 
#' \code{BDSimGeneral} will run until the number of total species generated, or, 
#' if \code{extOnly = TRUE}, the number of extant species at the end of the 
#' simulation, lies within the interval.
#' 
#' @param extOnly a boolean indicating whether \code{nFinal} should be taken as
#' the number of total or extant species during the simulation. If \code{TRUE},
#' \code{BDSimConstant} or \code{BDSimGeneral} will run until the number of extant
#' species lies within the \code{nFinal} interval. If \code{FALSE}, as default, it 
#' will run until the total number of species generated lies within that interval.
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
#' @return the return list of either \code{BDSimConstant} or
#' \code{BDSimGeneral}, which have the same elements, as follows
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
#' # we will showcase here some of the possible scenarios for diversification,
#' # touching on all the kinds of rates
#' 
#' ###
#' # consider first the simplest regimen, constant speciation and extinction
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
#' # run the simulation, making sure we have more than 1 species in the end
#' sim <- BDSim(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   # full phylogeny
#'   phy <- MakePhylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # now let us complicate speciation more, maybe a linear function
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation rate
#' p <- function(t) {
#'   return(0.05 + 0.005*t)
#' }
#' 
#' # extinction rate
#' q <- 0.08
#' 
#' # run the simulation, making sure we have more than 1 species in the end
#' sim <- BDSim(n0, p, q, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   # full phylogeny
#'   phy <- MakePhylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' # what if we want q to be a step function?
#' 
#' # list of extinction rates
#' qList <- c(0.08, 0.06, 0.07)
#' 
#' # list of shift times. Note qShifts could be c(40, 20, 5) for identical results
#' qShifts <- c(0, 20, 35)
#' 
#' # let us take a look at how MakeRate will make it a step function
#' q <- MakeRate(qList, fShifts = qShifts)
#' 
#' # and plot it
#' plot(seq(0, tMax, 0.1), q(seq(0, tMax, 0.1)), type = 'l',
#'      main = "Extintion rate as a step function", xlab = "Time (My)",
#'      ylab = "Rate (species/My)")
#' # note that this is slower than creating a step function with ifelse(), in this
#' # case q <- function(t) ifelse(t < 20, 0.08, ifelse(t < 35, 0.07, 0.08))
#' 
#' # also note that if done with ifelse(), the function must go from 0, instead of
#' # from tMax
#' 
#' # looking good, we will keep everything else the same
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # speciation
#' p <- function(t) {
#'   return(0.02 + 0.005*t)
#' }
#' 
#' # run the simulation. We can pass the step function directly, or just give
#' # a list of q and a list of shifts
#' sim <- BDSim(n0, p, qList, tMax, qShifts = qShifts, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # we can also supply a shape parameter to try age-dependent rates
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation - here note it is a Weibull scale
#' p <- 10
#' 
#' # speciation shape
#' pShape <- 2
#' 
#' # extinction
#' q <- 0.08
#' 
#' # run the simulation
#' sim <- BDSim(n0, p, q, tMax, pShape = pShape, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- MakePhylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # finally, we can also have a rate dependent on an environmental variable,
#' # like temperature data in RPANDA
#' 
#' if (requireNamespace("RPANDA", quietly = TRUE)) {
#'   
#'   # get temperature data from RPANDA
#'   data(InfTemp, package = "RPANDA")
#'   
#'   # initial number of species
#'   n0 <- 1
#'   
#'   # maximum simulation time
#'   tMax <- 40
#'   
#'   # speciation - a scale
#'   p <- function(t) {
#'     return(0.5 + 0.25*t)
#'   }
#'   
#'   # note the scale for the age-dependency can be a time-varying function
#'   
#'   # speciation shape
#'   pShape <- 1.5
#'   
#'   # extinction, dependent on temperature exponentially
#'   q <- function(t, env) {
#'     return(0.2*exp(0.01*env))
#'   }
#'   
#'   # need a variable to tell BDSim the extinction is environmentally dependent
#'   envQQ <- InfTemp
#'   
#'   # run the simulation
#'   sim <- BDSim(n0, p, q, tMax, pShape = pShape, envQQ = InfTemp, 
#'                nFinal = c(2, Inf))
#'   
#'   # we can plot the phylogeny to take a look
#'   if (requireNamespace("ape", quietly = TRUE)) {
#'     phy <- MakePhylo(sim)
#'     ape::plot.phylo(phy)
#'   }
#'   
#'   ###
#'   # one can mix and match all of these scenarios as they wish - age-dependency
#'   # and constant rates, age-dependent and temperature-dependent rates, etc. The
#'   # only combination that is not allowed is a vector rate and environmental
#'   # data, but one can get around that as follows
#'   
#'   # initial number of species
#'   n0 <- 1
#'   
#'   # speciation - a step function of temperature built using ifelse()
#'   p <- function(t, env) {
#'     ifelse(t < 20, env,
#'            ifelse(t < 30, env/2, 2*env/3))
#'   }
#'   
#'   # speciation shape
#'   pShape <- 2
#'   
#'   # environment variable to use - temperature
#'   envPP <- InfTemp
#'   
#'   # extinction
#'   q <- 0.15
#'   
#'   # maximum simulation time
#'   tMax <- 40
#'   
#'   # run the simulation
#'   sim <- BDSim(n0, p, q, tMax, pShape = pShape, envPP = envPP, 
#'                nFinal = c(2, Inf))
#'   
#'   # we can plot the phylogeny to take a look
#'   if (requireNamespace("ape", quietly = TRUE)) {
#'     phy <- MakePhylo(sim)
#'     ape::plot.phylo(phy)
#'   }
#' }
#' 
#' @name BDSim
#' @rdname BDSim
#' @export

BDSim <- function(n0, pp, qq, tMax, 
                  pShape = NULL, qShape = NULL, 
                  envPP = NULL, envQQ = NULL, 
                  pShifts = NULL, qShifts = NULL, 
                  nFinal = c(0, Inf), extOnly = FALSE,
                  fast = TRUE, trueExt=FALSE) {
  
  # if we have ONLY numbers for pp and qq, it is constant
  if ((is.numeric(pp) & length(pp) == 1) &
      (is.numeric(qq) & length(qq) == 1) &
       (is.null(c(pShape, qShape, envPP, envQQ, pShifts, qShifts)))) {
    p <- pp
    q <- qq
    
    # call BDSimConstant
    return(BDSimConstant(n0, p, q, tMax, nFinal, extOnly))
  }

  # else it is not constant
  # note even if pp or qq is constant this may call BDSimGeneral, since we
  # might have a shape parameter
  else {
    # use MakeRate to create the rates we want
    p <- MakeRate(pp, tMax, envPP, pShifts)
    q <- MakeRate(qq, tMax, envQQ, qShifts)

    # call BDSimGeneral
    return(BDSimGeneral(n0, p, q, tMax, pShape, qShape, 
                        nFinal, extOnly, fast, trueExt))
  }
}
