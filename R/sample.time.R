#' Constant and time-dependent rate species sampling
#' 
#' Generates a vector of occurrence times for species in a simulation using a
#' Poisson process. Allows for the Poisson rate to be (1) a constant or (2) a 
#' function of time. For fossil sampling dependent on species age in addition 
#' to absolute time, see \code{sample.age}. For a more general function that
#' considers these and other fossil sampling rate cases, see 
#' \code{sample.clade}.
#' 
#' Note that while the Poisson process occurs in forward time, we return (both 
#' in birth-death functions and here) results in backwards time, so that time 
#' is inverted using \code{tMax} both at the beginning and end of 
#' \code{sample.time}.
#'
#' @param S A vector representing the species numbers to be sampled. If 
#' \code{NULL}, as default, consider all species in the simulation.
#' 
#' @param rho Sampling rate (events per species per million years) over time. 
#' It can be a \code{numeric} describing a constant rate or a 
#' \code{function(t)} describing the variation in sampling over time. For more
#' flexibility on sampling, see \code{make.rate} for creating more complex 
#' rates. Note that \code{rho} should always be greater than or equal to zero.
#'
#' @inheritParams sample.clade
#'
#' @return A vector of occurrence times for that species.
#'
#' @author Bruno do Rosario Petrucci and Matheus Januario.
#'
#' @examples
#'
#' ###
#' # let us start with constant fossil sampling rate
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # sampling
#' rho <- 2
#' 
#' # time
#' time <-  seq(0, 10, by = 0.1)
#' 
#' # sample just the first species
#' occs <- sample.time(sim = sim, rho = rho, tMax = 10, S = 1)[[1]]
#' 
#' ###
#' # now let us try a step function
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#'
#' # we can create the sampling rate here from a few vectors
#' 
#' # rates
#' rList <-  c(0.5, 1, 1.5)
#' 
#' # rate shift times -  this could be c(10, 6, 2)
#' # and would produce the same function
#' rShifts <- c(0, 4, 8)
#' 
#' # create the rate to visualize it
#' rho <- make.rate(rList, tMax = 10, rateShifts = rShifts)
#' 
#' # time
#' time <-  seq(0, 10, by = 0.1)
#' 
#' # visualizing the plot from past to present
#' plot(x = time, y = rev(rho(time)), main = "Simulated preservation", 
#'      type = "l", xlab = "Time (Mya)", ylab = "Rate (events/species/My)",
#'      xlim = c(10, 0))
#' 
#' \dontrun{
#' # sample first two species
#' occs <- sample.time(sim = sim, rho = rho, tMax = 10, S = 1:2)
#' }
#' 
#' ###
#' # we can create a step function using ifelse as well
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # preservation function
#' rho <- function(t) {
#'   ifelse(t < 4, 0.5,
#'          ifelse(t < 8, 1, 1.5))
#' }
#' # note how this function should be exactly the same as the previous one
#' 
#' # time
#' time <-  seq(0, 10, by = 0.1)
#' 
#' # visualizing the plot from past to present
#' plot(x = time, y = rev(rho(time)), main = "Simulated preservation", 
#'      type = "l", xlab = "Time (Mya)", ylab = "Rate (events/species/My)",
#'      xlim = c(10, 0))
#' 
#' # sample all species
#' occs <- sample.time(sim = sim, rho = rho, tMax = 10)
#' 
#' ###
#' # finally, we could generate sampling dependent on a time-series
#' # in this case, an environmental variable, temperature
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # get temperature data
#' data(temp)
#' 
#' # preservation function dependent on temperature
#' r_t <-  function(t, temp) {
#'   return(0.1*temp)
#' }
#' 
#' # final preservation
#' rho <- make.rate(r_t, tMax = tMax, envRate = temp)
#' 
#' # visualizing the plot from past to present
#' plot(x = time, y = rev(rho(time)), main = "Simulated preservation", 
#'      type = "l", xlab = "Time (Mya)", ylab = "Rate (events/species/My)",
#'      xlim = c(10, 0))
#' 
#' # sample the first species
#' occs <- sample.time(sim = sim, rho = rho, tMax = 10, S = 1)
#' 
#' # after presenting the possible models, we can consider how to
#' # create mixed models, where the dependency changes over time
#' 
#' ###
#' # consider sampling that becomes temperature dependent
#' # in the middle of the simulation
#' 
#' # set seed
#' set.seed(1)
#'
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#'
#' # get the temperature data
#' data(temp)
#'
#' # preservation function dependent on t and temperature
#' r_t <-  function(t, temp) {
#'   return(
#'     ifelse(t < 5, 1 - 0.25*t,
#'            0.1*temp)
#'   )
#' }
#' 
#' # final preservation
#' rho <- make.rate(r_t, tMax = tMax, envRate = temp)
#' 
#' # visualizing the plot from past to present
#' plot(x = time, y = rev(rho(time)), main = "Simulated preservation", 
#'      type = "l", xlab = "Time (Mya)", ylab = "Rate (events/species/My)",
#'      xlim = c(10, 0))
#' 
#' # sample
#' occs <- sample.time(sim = sim, rho = rho, tMax = 10, S = 1)
#' 
#' ###
#' # we can also change the environmental variable
#' # halfway into the simulation
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # get the temperature data
#' data(temp)
#' 
#' # temperature-dependent preservation
#' r_t1 <- function(t, temp) {
#'   return(0.08*temp)
#' }
#' 
#' # make first function
#' rho1 <- make.rate(r_t1, tMax = tMax, envRate = temp)
#' 
#' # get the co2 data
#' data(co2)
#' 
#' # co2-dependent preservation
#' r_t2 <- function(t, co2) {
#'   return(1.5 - 0.05*co2)
#' }
#' 
#' # make second function
#' rho2 <- make.rate(r_t2, tMax = tMax, envRate = co2)
#' 
#' # final preservation function
#' rho <- function(t) {
#'   ifelse(t < 5, rho1(t), rho2(t))
#' }
#' 
#' # visualizing the plot from past to present
#' plot(x = time, y = rev(rho(time)), main = "Simulated preservation", 
#'      type = "l", xlab = "Time (Mya)", ylab = "Rate (events/species/My)",
#'      xlim = c(10, 0))
#' 
#' # sample
#' occs <- sample.time(sim = sim, rho = rho, tMax = 10, S = 1)
#' 
#' @name sample.time
#' @rdname sample.time
#' @export

sample.time <- function(sim, rho, tMax, S = NULL) {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # if S is NULL, make it all species in the simulation
  if (is.null(S)) {
    S <- 1:length(sim$TE)
  }
  
  # if there are species in S not in the simulation, error
  if (!(all(S %in% 1:length(sim$TE)))) {
    stop("S must contain only integers between 1 and the total 
         number of species in sim")
  }
  
  # get speciation and extinction times
  TE <- sim$TE
  TS <- sim$TS
  
  # make TE sensible
  TE[sim$EXTANT] <- 0
  
  # invert them since simulation goes from 0 to tMax
  TE <- tMax - TE
  TS <- tMax - TS
  
  # create return
  res <- list()
  
  # for each species
  for (s in S) {
    # start when the species was born, end when it died
    now <- TS[s]
    End <- TE[s]
    
    # initialize vector
    sampled <- c()
    
    # while we have not reached the time of extinction
    while (now < End) {
      # take the waiting time for sampling, using rexp.var()
      WaitTimeR <- ifelse(is.numeric(rho), 
                          ifelse(rho > 0, rexp(1, rho), Inf),
                          ifelse(rho(now) > 0, rexp.var(1, rho, now, tMax,
                                                        fast = TRUE), Inf))
      
      # advance in time
      now <- now + WaitTimeR
      
      # if sampling comes after extinction, we don't include this occurrence
      if (now > End) break
      
      # append to the vector
      sampled <- c(sampled, now)
    }
    
    # finally, we invert time so that it goes from tMax to 0
    sampled <- tMax - sampled
    
    # append to result
    res[[s]] <- sampled
  }
  
  return(res)
}
