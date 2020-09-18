#' Constant and time-dependent rate species sampling
#' 
#' Generates a vector of occurrence times for a species in a simulation using a
#' a Poisson process. Allows for the Poisson rate to be (1) a constant or (2) a 
#' function of time. For sampling of more than one species and/or taking into 
#' account species age in addition to absolute time, see \code{sample.clade} and
#' \code{sample.general}. \code{sample.clade} also allows for more flexibility
#' options, see \code{make.rate}.
#' Note that while the Poisson process occurs in forward time, we return (both in
#' birth-death functions and here) results in backwards time, so that time is
#' inverted using \code{tMax} both at the beginning and end of 
#' \code{sample.species}.
#'
#' @param S The species number to be sampled. Since \code{sample.species} will be 
#' called by a wrapper using \code{lapply}, it is through \code{S} that we apply
#' this function.
#' 
#' @param rho Sampling rate (per species per million years) over time. It can be
#' a \code{numeric} describing a constant rate or a \code{function(t)} describing
#' the variation in sampling over time. For more flexibility on sampling, see
#' \code{make.rate} for creating more complex rates. Note that \code{rho} should
#' always be greater than or equal to zero.
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
#' # let us start with a linear increase in preservation rate
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation was short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # preservation function
#' rho <- function(t) {
#'   return(1 + 0.25*t)
#' }
#' 
#' # time
#' t <- seq(0, 10, by = 0.1)
#' 
#' # visualizing from the past to the present
#' plot(x = t, y = rev(rho(t)), main="Simulated preservation", type = "l",
#'      xlab = "Mya", ylab = "preservation rate",
#'      xlim = c(10, sim$TE[1]))
#' 
#' # sample
#' occs <- sample.species(sim = sim, rho = rho, tMax = 10, S = 1)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, sim$TE[1]),
#'      xlab = "Mya")
#' lines(t, rev(rho(t)))
#' 
#' ###
#' # now let us try a step function
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation was short lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # we can create the sampling rate here from a few vectors
#' 
#' # rates
#' rList <- c(1, 3, 0.5)
#' 
#' # rate shift times -  this could be c(10, 6, 2)
#' # and would produce the same function
#' rShifts <- c(0, 4, 8)
#' 
#' # create the rate to visualize it
#' r <- make.rate(rList, tMax = 10, rateShifts = rShifts)
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
#' occs <- sample.species(sim = sim, rho = r, tMax = 10, S = 1)
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
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation was short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # preservation function
#' rho <- function(t) {
#'   ifelse(t < 4, 1,
#'          ifelse(t < 8, 3, 0.5))
#' }
#' # note how this function should be exactly the same as the previous one
#' 
#' # time
#' t <- seq(0, 10, by = 0.1)
#' 
#' # visualizing the plot from past to present
#' plot(x = t, y = rev(rho(t)), main = "Simulated preservation", type = "l",
#'      xlab = "Mya", ylab = "preservation rate",
#'      xlim = c(10, sim$TE[1]))
#' 
#' # sample
#' occs <- sample.species(sim = sim, rho = rho, tMax = 10, S = 1)
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
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation was short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # preservation function dependent on temperature
#' r_t <- function(t, env) {
#'   return(0.25*env)
#' }
#' 
#' # get the temperature data
#' data(temp)
#' 
#' # final preservation
#' rho <- make.rate(r_t, envRate = temp)
#' 
#' # visualizing the plot from past to present
#' plot(x = t, y = rev(r(t)), main = "Simulated preservation", type = "l",
#'      xlab = "Mya", ylab = "preservation rate",
#'      xlim = c(10, ifelse(is.na(sim$TE[1]), 0, sim$TE[1])))
#' 
#' # sample
#' occs <- sample.species(sim = sim, rho = rho, tMax = 10, S = 1)
#' 
#' # check histogram
#' hist(occs,
#'      xlim = c(10, ifelse(is.na(sim$TE[1]), 0, sim$TE[1])),
#'      xlab = "Mya")
#' lines(t, rev(rho(t)))
#' 
#' @name sample.species
#' @rdname sample.species
#' @export

sample.species <- function(sim, rho, tMax, S) {
  # get speciation and extinction times
  TE <- sim$TE
  TS <- sim$TS
  
  # make TE sensible
  TE[sim$EXTANT] <- 0
  
  # invert them since simulation goes from 0 to tMax
  TE <- tMax - TE
  TS <- tMax - TS
  
  # start when the species was born, end when it died
  now <- TS[S]
  End <- TE[S]

  # make rho a function if it isn't
  if (is.numeric(rho)) {
    r <- rho
    rho <- Vectorize(function(t) r)
  }
  
  # initialize vector
  sampled <- c()

  # while we have not reached the time of extinction
  while (now < End) {
    # take the waiting time for sampling, using rexp.var()
    WaitTimeR <- ifelse(rho(now) > 0, rexp.var(1, rho, now, tMax), Inf)

    # advance in time
    now <- now + WaitTimeR

    # if sampling comes after extinction, we don't include this occurrence
    if (now > End) break

    # append to the vector
    sampled <- c(sampled, now)
  }

  # finally, we invert time so that it goes from tMax to 0
  sampled <- tMax - sampled

  return(sampled)
}
