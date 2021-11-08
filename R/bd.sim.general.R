#' Non-constant rate Birth-Death simulation
#'
#' Simulates a species birth-death process with general rates for any number of
#' starting species. Allows for the speciation/extinction rate to be (1) a 
#' constant, or (2) a function of time. Allows for constraining results on the 
#' number of species at the end of the simulation, either total or extant. The 
#' function can also take an optional shape argument to generate age-dependence
#' on speciation and/or extinction, assuming a Weibull distribution as a model 
#' of age-dependence. Returns a \code{sim} object (see \code{?sim}). It may 
#' return calculated extinction times or simply information on whether species 
#' lived after the maximum simulation time, depending on input. For constant 
#' rate simulations, see \code{bd.sim.constant}. For a function that unites all
#' scenarios, see \code{bd.sim}. \code{bd.sim} also allows for extra inputs, 
#' creating a time-dependent only rate internally through \code{make.rate}. For
#' similar flexibility, use \code{make.rate} to generate the desired rate.
#' Please note while time runs from \code{0} to \code{tMax} in the simulation,
#' it returns speciation/extinction times as \code{tMax} (origin of the group) 
#' to \code{0} (the "present" and end of simulation), so as to conform to other
#' packages in the literature.
#'
#' @inheritParams bd.sim
#' 
#' @param lambda Function that defines the speciation rate over time. It will 
#' either be interpreted the rate of an exponential distribution, or the scale 
#' of a Weibull distribution if \code{lShape != NULL}. Can be constant, to 
#' allow for mixing of constant and non-constant rates. One can use constructs 
#' such as \code{ifelse()} to create rates whose underlying model change over 
#' time (see the last examples). Note that \code{lambda} should always be 
#' greater than or equal to zero.
#'
#' @param mu Similar to above, but for the extinction rate.
#' 
#' Note: rates should be considered as running from \code{0} to \code{tMax}, as
#' the simulation runs in that direction even though the function inverts 
#' speciation and extinction times before returning.
#' 
#' Note: this function is meant to be called by \code{bd.sim}, so it neither
#' allows for as much flexibility, nor calls \code{make.rate}. If the user 
#' wishes to use \code{bd.sim.general} with environmental or step-function 
#' rates, they can generate the rate with \code{make.rate} and supply it to the
#' function.
#'
#' @return A \code{sim} object, containing extinction times, speciation times,
#' parent, and status information for each species in the simulation. See 
#' \code{?sim}.
#'
#' @author Bruno do Rosario Petrucci.
#'
#' @examples
#'
#' # we will showcase here some of the possible scenarios for diversification,
#' # touching on all the kinds of rates
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
#' lambda <- 0.11
#' 
#' # extinction
#' mu <- 0.08
#' 
#' # set a seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # we can complicate things further with a linear function as a rate
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # make a vector for time
#' time <- seq(0, tMax, 0.1)
#' 
#' # speciation rate
#' lambda <- function(t) {
#'   return(0.05 + 0.005*t)
#' }
#' 
#' # extinction rate
#' mu <- 0.1
#' 
#' # set a seed
#' set.seed(4)
#' 
#' # run the simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, nFinal = c(2, Inf))
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
#' # speciation
#' lambda <- function(t) {
#'   return(0.02 + 0.005*t)
#' }
#' 
#' # vector of extinction rates
#' mList <- c(0.09, 0.08, 0.1)
#' 
#' # vector of shift times. Note mShifts could be c(40, 20, 5) for identical 
#' # results
#' mShifts <- c(0, 20, 35)
#' 
#' # let us take a look at how make.rate will make it a step function
#' mu <- make.rate(mList, tMax = tMax, rateShifts = mShifts)
#' 
#' # and plot it
#' plot(seq(0, tMax, 0.1), rev(mu(seq(0, tMax, 0.1))), type = 'l',
#'      main = "Extintion rate as a step function", xlab = "Time (My)",
#'      ylab = "Rate (species/My)", xlim = c(tMax, 0))
#' 
#' # looking good, we will keep everything else the same
#' 
#' # a different way to define the same extinction function
#' mu <- function(t) {
#'   ifelse(t < 20, 0.09,
#'          ifelse(t < 35, 0.08, 0.1))
#' }
#' 
#' # set seed
#' set.seed(2)
#' 
#' # run the simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, nFinal = c(2, Inf))
#' # we could instead have used mu made with make.rates
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # another feature to add is age dependency
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' lambda <- 0.1
#' 
#' # extinction - note that since this is a Weibull scale,
#' # the unites are my/events/lineage, not events/lineage/my
#' mu <- 10
#' 
#' # extinction shape
#' mShape <- 2
#' 
#' # set seed
#' set.seed(3)
#' 
#' # run simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, mShape = mShape, 
#'                       nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # scale can be a time-varying function  
#'  
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' lambda <- 0.15
#' 
#' # extinction - note that since this is a Weibull scale,
#' # the unites are my/events/lineage, not events/lineage/my
#' mu <- function(t) {
#'   return(8 + 0.05*t)
#' }
#' 
#' # extinction shape
#' mShape <- 2
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, mShape = mShape, 
#'                       nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # and shape can also vary with time
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' lambda <- 0.15
#' 
#' # extinction - note that since this is a Weibull scale,
#' # the unites are my/events/lineage, not events/lineage/my
#' mu <- 5
#' 
#' # extinction shape
#' mShape <- function(t) {
#'   return(8 + 0.05*t)
#' }
#' 
#' # set seed
#' set.seed(2)
#' 
#' # run simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, mShape = mShape, 
#'                       nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # finally, we could have environmental dependency on a rate
#' 
#' # get the temperature data - see ?temp for information on the data set
#' data(temp)
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # temperature-dependent speciation
#' l_t <- function(t, temp) {
#'  return(0.025*exp(0.1*temp))
#' }
#' 
#' # need a data frame describing the temperature at different times
#' envL <- temp
#' 
#' # speciation - note that bd.sim does this internally when supplied
#' # with the equivalent of l_t and envL
#' lambda <- make.rate(l_t, tMax = tMax, envRate = envL)
#' 
#' # extinction
#' mu <- 0.08
#' 
#' # set seed
#' set.seed(3)
#'
#' # run simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # one can mix and match all of these scenarios as they wish - age-dependency
#' # and constant rates, age-dependent and temperature-dependent rates, etc. 
#' # the only combination that is not allowed by make.rate is a step function 
#' # rate and environmental data, but one can get around that as follows
#' 
#' # get the temperature data - see ?temp for information on the data set
#' data(temp)
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # speciation - a step function of temperature built using ifelse()
#' # note that this creates two shifts for lambda, for a total of 3 values
#' # throughout the simulation
#' l_t <- function(t, env) {
#'   ifelse(t < 20, env,
#'          ifelse(t < 30, env / 4, env / 3))
#' }
#' 
#' # speciation shape
#' lShape <- 2
#' 
#' # environment variable to use - temperature
#' envL <- temp
#' 
#' # speciation
#' lambda <- make.rate(l_t, tMax = tMax, envRate = envL)
#' 
#' # plot it to take a look
#' plot(seq(0, tMax, 0.1), rev(lambda(seq(0, tMax, 0.1))),
#'      main = "Speciation scale varying with temperature", xlab = "Time (My)",
#'      ylab = "Scale", type = 'l', xlim = c(tMax, 0))
#' 
#' # extinction
#' mu <- 0.1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # set seed
#' set.seed(1)
#' 
#' \dontrun{
#' # run the simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, lShape = lShape,
#'                       nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' }
#' 
#' # after presenting the possible models, we can consider how to
#' # create mixed models, where the dependency changes over time
#' 
#' ###
#' # consider speciation that becomes environment dependent
#' # in the middle of the simulation
#' 
#' # get the temperature data - see ?temp for information on the data set
#' data(temp)
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # time and temperature-dependent speciation
#' l_t <- function(t, temp) {
#'   return(
#'     ifelse(t < 20, 0.1 - 0.005*t,
#'            0.05 + 0.1*exp(0.02*temp))
#'   )
#' }
#' 
#' # speciation
#' lambda <- make.rate(l_t, tMax = tMax, envRate = temp)
#' 
#' # extinction
#' mu <- 0.1
#' 
#' # set seed
#' set.seed(4)
#' 
#' \dontrun{
#' # run simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' }
#' 
#' ###
#' # we can also change the environmental variable
#' # halfway into the simulation
#' 
#' # get the temperature data - see ?temp for information on the data set
#' data(temp)
#' 
#' # same for co2 data (and ?co2)
#' data(co2)
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' lambda <- 0.1
#' 
#' # temperature-dependent extinction
#' m_t1 <- function(t, temp) {
#'   return(0.05 + 0.1*exp(0.02*temp))
#' }
#' 
#' # make first function
#' mu1 <- make.rate(m_t1, tMax = tMax, envRate = temp) 
#' 
#' # co2-dependent extinction
#' m_t2 <- function(t, co2) {
#'   return(0.02 + 0.14*exp(0.01*co2))
#' }
#' 
#' # make second function
#' mu2 <- make.rate(m_t2, tMax = tMax, envRate = co2)
#' 
#' # final extinction function
#' mu <- function(t) {
#'   ifelse(t < 20, mu1(t), mu2(t))
#' }
#' 
#' # set seed
#' set.seed(3)
#' 
#' # run simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' # note one can also use this mu1 mu2 workflow to create a rate
#' # dependent on more than one environmental variable, by decoupling
#' # the dependence of each in a different function and putting those
#' # together
#' 
#' ###
#' # finally, one could create an extinction rate that turns age-dependent
#' # in the middle, by making shape time-dependent
#'
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' lambda <- 0.15
#' 
#' # extinction - note that since this is a Weibull scale,
#' # the unites are my/events/lineage, not events/lineage/my
#' mu <- function(t) {
#'   return(8 + 0.05*t)
#' }
#' 
#' # extinction shape
#' mShape <- function(t) {
#'   return(
#'     ifelse(t < 30, 1, 2)
#'   )
#' }
#' 
#' # set seed
#' set.seed(3)
#' 
#' # run simulation
#' sim <- bd.sim.general(n0, lambda, mu, tMax, mShape = mShape,
#'                       nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # note nFinal has to be sensible
#' \dontrun{
#' # this would return a warning, since it is virtually impossible to get 100
#' # species at a process with diversification rate -0.09 starting at n0 = 1
#' sim <- bd.sim.general(1, lambda = 0.01, mu = 1, tMax = 100, 
#'                       nFinal = c(100, Inf))
#' }
#'
#' @name bd.sim.general
#' @rdname bd.sim.general
#' @export
#' 

bd.sim.general <- function(n0, lambda, mu, tMax, 
                         lShape = NULL, mShape = NULL,
                         nFinal = c(0, Inf), nExtant = c(0, Inf),
                         trueExt = FALSE) {
  # error check - rate cannot be negative
  if ((is.numeric(lambda))) {
    if (lambda < 0) {
      stop("speciation rate cannot be negative")
    }
  }
  
  else {
    if (optimize(lambda, interval = c(0, 1e10))$objective < 0) {
      stop("speciation rate cannot be negative at any point in time")
    }
  }
  
  if ((is.numeric(mu))) {
    if (mu < 0) {
      stop("extinction rate cannot be negative")
    }
  }
  
  else {
    if (optimize(mu, interval = c(0, 1e10))$objective < 0) {
      stop("extinction rate cannot be negative at any point in time")
    }
  }
  
  # check that n0 is not negative
  if (n0 <= 0) {
    stop("initial number of species must be positive")
  }
  # check nFinal is sensible - two numbers, maximum >=1
  if ((length(nFinal) != 2) || (typeof(nFinal) != "double")) {
    stop("nFinal must be a vector with a minimum and maximum number 
         of species")
  } else if (max(nFinal) < 1) {
    stop("nFinal must have a maximum number of species greater than 0")
  } else {
    # if everything is good, make sure it's sorted
    nFinal <- sort(nFinal)
  }
  
  # similarly for nExtant
  if ((length(nExtant) != 2) || (typeof(nExtant) != "double")) {
    stop("nExtant must be a vector with a minimum and maximum number 
         of species")
  } else if (max(nExtant) < 0) {
    stop("nExtant must have a maximum number of species greater 
         than or equal to 0")
  } else {
    # if everything is good, make sure it's sorted
    nExtant <- sort(nExtant)
  }
  # if shape is not null, make scale a function to facilitate checking
  if (!is.null(lShape)) {
    message("since lShape is not null, lambda will be a Weibull scale and
            therefore correspond to 1/rate. See ?bd.sim or ?bd.sim.general")
    
    if (is.numeric(lambda)) {
      l <- lambda
      lambda <- Vectorize(function(t) l)
    }
    
    # check that it is never <= 0
    if (is.numeric(lShape)) {
      if (lShape <= 0) {
        stop("lShape may be nonpositive. It must always be >0")
      }
    }
    
    else {
      if (optimize(lShape, interval = c(0, 1e10))$objective < 0.01) {
        stop("lShape may be nonpositive. It must always be >0")
      }
    }
  }  
  
  if (!is.null(mShape)) {
    message("since mShape is not null, mu will be a Weibull scale and therefore
            correspond to 1/rate. See ?bd.sim or ?bd.sim.general")
    
    if (is.numeric(mu)) {
      m <- mu
      mu <- Vectorize(function(t) m)
    }
    
    # check that it is never <= 0
    if (is.numeric(mShape)) {
      if (mShape <= 0) {
        stop("mShape may be nonpositive. It must always be >0")
      }
    }
    
    else {
      if (optimize(mShape, interval = c(0, 1e10))$objective < 0.01) {
        stop("mShape may be nonpositive. It must always be >0")
      }
    }
  }
 
  # initialize test making sure while loop runs
  inBounds <- FALSE
  
  # counter to make sure the nFinal is achievable
  counter <- 1

  while (!inBounds) {
    # create vectors to hold times of speciation, extinction, 
    # parents and status
    TS <- rep(0, n0)
    TE <- rep(NA, n0)
    parent <- rep(NA, n0)
    isExtant <- rep(TRUE, n0)
  
    # initialize species count
    sCount <- 1
  
    # while we have species to be analyzed still
    while (length(TE) >= sCount) {
      # start at the time of speciation of sCount
      tNow <- TS[sCount]

      # find the waiting time using rexp.var if lambda is not constant
      waitTimeS <- ifelse(
        is.numeric(lambda), ifelse(lambda > 0, rexp(1, lambda), Inf),
        ifelse(lambda(tNow) > 0, 
               rexp.var(1, lambda, tNow, tMax, lShape, TS[sCount],
                        fast = TRUE), Inf))
      waitTimeE <- ifelse(
        is.numeric(mu), ifelse(mu > 0, rexp(1, mu), Inf),
        ifelse(mu(tNow) > 0,
               rexp.var(1, mu, tNow, tMax, mShape, TS[sCount],
                        fast = !trueExt), Inf))
      # fast is true since we do not record speciations after
      # the extinction anyway
  
      tExp <- tNow + waitTimeE
  
      # while there are fast enough speciations before the species 
      # goes extinct,
      while ((tNow + waitTimeS) < min(tExp, tMax)) {
  
        # advance to the time of speciation
        tNow <- tNow + waitTimeS
  
        # add new times to the vectors
        TS <- c(TS, tNow)
        TE <- c(TE, NA)
        parent <- c(parent, sCount)
        isExtant <- c(isExtant, TRUE)
  
        # get a new speciation waiting time, and include it in the vector
        waitTimeS <- ifelse(
          is.numeric(lambda), ifelse(lambda > 0, rexp(1, lambda), Inf),
          ifelse(lambda(tNow) > 0, 
                 rexp.var(1, lambda, tNow, tMax, lShape, TS[sCount],
                          fast = TRUE), Inf))
        # fast is true since we do not record speciations after
        # the extinction anyway
      }
  
      # reached the time of extinction
      tNow <- tExp
  
      # if trueExt is true or the species went extinct before tMax,
      # record it. If both are false record it as NA
      TE[sCount] <- ifelse(tNow < tMax | trueExt, tNow, NA)
      
      # record the extinction -
      isExtant[sCount] <- ifelse(is.na(TE[sCount]) | TE[sCount] > tMax,
                                 TRUE, FALSE)
  
      # next species
      sCount <- sCount + 1
    }
  
    # now we invert TE and TS so time goes from tMax to 0
    TE <- tMax - TE
    TS <- tMax - TS

    # check whether we are in bounds
    inBounds <- (length(TE) >= nFinal[1]) &&
      (length(TE) <= nFinal[2]) &&
      (sum(isExtant) >= nExtant[1]) &&
      (sum(isExtant) <= nExtant[2])
    
    # if we have ran for too long, stop
    counter <- counter + 1
    if (counter > 100000) {
      warning("This value of nFinal took more than 100000 simulations 
              to achieve")
      return(NA)
    }
  }
  
  # create the return
  sim <- list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant)
  class(sim) <- "sim"
  
  return(sim)
}
