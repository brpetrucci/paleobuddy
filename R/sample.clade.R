#' General rate species sampling
#' 
#' Generates occurrence times or time ranges (as most empirical fossil 
#' occurrences) for each of the desired species using a Poisson process. Allows
#' for the Poisson rate to be (1) a constant, (2) a function of time, (3) a 
#' function of time and a time-series (usually environmental) variable, or (4) 
#' a vector of numbers (rates in a step function). Allows for age-dependent 
#' sampling with a parameter for a distribution representing the expected 
#' occurrence number over a species duration. Allows for further flexibility in
#' rates by a shift times vector and environmental matrix parameters. 
#' Optionally takes a vector of time bins representing geologic periods, if the
#' user wishes occurrence times to be represented as a range instead of true 
#' points. See \code{sample.time} - absolute time-dependent sampling only - 
#' and \code{sample.age} - time and/or age-dependent sampling - for more 
#' information.
#'
#' @param bins A vector of time intervals corresponding to geological time ranges.
#' If it is not supplied, \code{seq(tMax, 0, -0.1)} is used.
#' 
#' @inheritParams sample.age
#'
#' @param rho Sampling rate (per species per million years) over time. It can be 
#' a \code{numeric} describing a constant rate, a \code{function(t)} describing 
#' the variation in sampling over time \code{t}, a \code{function(t, env)} 
#' describing the variation in sampling over time following both time AND 
#' an environmental variable (please see \code{envR} for details), or a 
#' \code{vector} containing rates that correspond to each rate between sampling
#' rate shift times times (please see \code{rShifts}). Note that \code{rho} should
#' should always be greater than or equal to zero.
#' 
#' @param envR A data frame containing time points and values of an environmental
#' variable, like temperature, for each time point. This will be used to create
#' a sampling rate, so \code{rho} must be a function of time and said variable
#' if \code{envR} is not \code{NULL}. Note \code{paleobuddy} has two 
#' environmental data frames, \code{temp} and \code{co2}. See \code{RPANDA} for
#' more examples.
#'
#' @param rShifts Vector of rate shifts. First element must be the starting
#' time for the simulation (\code{0} or \code{tMax}). It must have the same length
#' as \code{lambda}. \code{c(0, x, tMax)} is equivalent to 
#' \code{c(tMax, tMax - x, 0)} for the purposes of \code{make.rate}.
#'
#' @param returnTrue If set to \code{FALSE}, it will contain the occurrence
#' times as ranges. In this way, we simulate the granularity presented by
#' empirical fossil records. If \code{returnTrue} is \code{TRUE}, this is ignored.
#' 
#' @param returnAll If set to \code{TRUE}, returns both the true sampling time and
#' age ranges. Default is \code{FALSE}
#'
#' @return A \code{data.frame} containing species names/numbers, whether each 
#' species is extant, and the true occurrence times of each fossil, a range of 
#' occurrence times based on \code{bins}, or both.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @examples
#' 
#' ###
#' # we can start with a constant case
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # sampling rate
#' rho <- 2
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' # find the occurrence data frame
#' dt <- sample.clade(sim, rho, tMax = 10, bins = bins, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[i]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)",
#'        xlim = c(sim$TS[i], sim$TE[i]))
#' }
#' 
#' ###
#' # sampling can be any function of time 
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(3 - 0.15*t)
#' }
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' # find the occurrence data frame
#' dt <- sample.clade(sim, rho, tMax = 10, bins = bins, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[i]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)",
#'        xlim = c(sim$TS[i], sim$TE[i]))
#' }
#' 
#' ###
#' # now we can try a step function rate
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # we will use the less efficient method of creating a step function
#' # one could instead use ifelse()
#' 
#' # rates vector
#' rList <- c(1, 2, 0.5)
#' 
#' # rate shifts vector
#' rShifts <- c(0, 4, 8)
#' 
#' # make it a function so we can plot it
#' rho <- make.rate(rList, 10, rateShifts = rShifts)
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' # find the occurrence data frame
#' dt <- sample.clade(sim, rList, rShifts = rShifts, 
#'                    tMax = 10, bins = bins, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[i]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)",
#'        xlim = c(sim$TS[i], sim$TE[i]))
#' }
#' 
#' ###
#' # finally, sample.clade also accepts an environmental variable
#' 
#' # get temperature data
#' data(temp)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # make temperature the environmental dependency of r
#' envR <- temp
#' 
#' # we can then make sampling dependent on the temperature
#' r_t <- function(t, env) {
#'   return(0.5*env)
#' }
#' 
#' # make it a function so we can plot it
#' rho <- make.rate(r_t, tMax = tMax, envRate = envR)
#' 
#' # let us check that rho is high enough to see a pattern
#' plot(1:10, rho(1:10), type = 'l', main = "Sampling rate",
#'      xlab = "My", ylab = "rho")
#' 
#' # the resolution of the fossil data set:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' # find the occurrence data frame
#' dt <- sample.clade(sim, r_t, tMax = 10, envR = envR,
#'                    bins = bins, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#' 
#'   # check the histogram
#'   hist(mids[dt$Species == ids[i]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)",
#'        xlim = c(sim$TS[i], sim$TE[i]))
#' }
#' 
#' # we will now do some examples with age-dependent rates. For more details,
#' # check sample.age.
#' 
#' ###
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # here we will use the PERT function. It is described in:
#' # Silvestro et al 2014
#' 
#' # preservation function
#' dPERT <- function(t, s, e, sp, a = 3, b = 3, log = FALSE) {
#'   
#'   # check if it is a valid PERT
#'   if (e >= s) {
#'     message("There is no PERT with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # find the valid and invalid times
#'   id1 <- which(t <= e | t >= s)
#'   id2 <- which(!(t <= e | t >= s))
#'   t <- t[id2]
#'   
#'   # initialize result vector
#'   res <- vector()
#'   
#'   # if user wants a log function
#'   if (log) {
#'     # invalid times get -Inf
#'     res[id1] <- -Inf
#'     
#'     # valid times calculated with log
#'     res[id2] <- log(((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b)))
#'   }
#'   # otherwise
#'   else{
#'     res[id1] <- 0
#'     
#'     res[id2] <- ((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b))
#'   }
#'   
#'   return(res)
#' }
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' dt <- sample.clade(sim, rho = 3, tMax = 10, bins = bins,
#'                    adFun = dPERT, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[[i]]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)", probability = TRUE)
#'   
#'   # expected curve
#'   curve(dPERT(x, s = sim$TS[sp], e = sim$TE[sp], sp = sp), from = sim$TE[sp],
#'         to = sim$TS[sp], add = TRUE, col = "red", n = 100)
#' }
#' # we provide curves for comparison here, but remember the low sample sizes 
#' # (and bins) may affect the quality of the fit
#' 
#' ###
#' # now, a hat-shaped increase through the duration of a species dependent on two
#' # parameters
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # get parameters 
#' 
#' # a random point inside each lineage's duration
#' par <- runif (n = length(sim$TE), min = sim$TE, max = sim$TS)
#' 
#' # a distance between "par" and the lineage's duration middle
#' par1 <- (((sim$TS - sim$TE) / 2) + sim$TE) - par
#' 
#' # preservation function in respect to age, with the "mode" of the triangle
#' # being dependent on par and par1
#' dTRImod <- function(t, s, e, sp) {
#'   
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # md depends on the two parameters
#'   md <- par[sp] + par1[sp]
#'   
#'   # check that md is valid
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t > md & t <= s)
#'   id4 <- which(!(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   res <- vector()
#'   
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e)*(md - e))
#'   res[id2] <- 2 / (s - e)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e)*(s - md))
#'   res[id4] <- 0
#'   
#'   return(res)
#' }
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' dt <- sample.clade(sim, rho = 4, tMax = 10, bins = bins,
#'                    adFun = dTRImod, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[[i]]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)", probability = TRUE)
#'   
#'   # expected curve
#'   curve(dTRImod(x, e = sim$TE[sp], s = sim$TS[sp], sp = sp),from = sim$TE[sp],
#'         to = sim$TS[sp], add = TRUE, col="red", n = 100)
#' }
#' # we provide curves for comparison here, but remember the low sample sizes 
#' # (and bins) may affect the quality of the fit
#' 
#' ###
#' # let us keep everything from the last example, but
#' # having a time-dependent sampling rate rho
#' 
#' # in this case, the function finds the number of
#' # occurrences using rho, and their distribution
#' # using a normalized version of rho * adFun
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(2 + 0.1*t)
#' }
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.2, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # get parameters
#' 
#' # a random point inside each lineage's duration
#' par <- runif (n = length(sim$TE), min = sim$TE, max = sim$TS)
#' 
#' # a distance between "par" and the lineage's duration middle
#' par1 <- (((sim$TS - sim$TE) / 2) + sim$TE) - par
#' 
#' # preservation function in respect to age, with the "mode" of the triangle
#' # being exactly at the last quarter of the duration of EACH lineage.
#' dTRImod <- function(t, s, e, sp) {
#'   
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # md depends on the two parameters
#'   md <- par[sp] + par1[sp]
#'   
#'   # check that md is valid
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t > md & t <= s)
#'   id4 <- which(!(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   res <- vector()
#'   
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e)*(md - e))
#'   res[id2] <- 2 / (s - e)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e)*(s - md))
#'   res[id4] <- 0
#'   
#'   return(res)
#' }
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' dt <- sample.clade(sim, rho = rho, tMax = 10, bins = bins,
#'                    adFun = dTRImod, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[[i]]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)", probability = TRUE)
#'   
#'   # expected curve
#'   tt <- seq(from = sim$TE[sp], to = sim$TS[sp], by = 0.01)
#'   lines(x = tt, y = dTRImod(tt, s = sim$TS[sp], e = sim$TE[sp], sp = sp))
#'   
#'   # need tMax to invert rho
#'   tMax <- 10
#'   
#'   # getting expected values from age + time dependences:
#'   Pres_time_adpp <- function(t, s, e, sp, ...) {
#'     # correction of scale
#'     rhoMod <- function(t) {
#'       return(rho(tMax - t))
#'     }
#'     return(rhoMod(t)*dTRImod(t = t, s = s, e = e, sp = sp, ...))
#'   }
#'   
#'   # normalizing
#'   Pres_time_adppNorm <- function(t, s, e, sp, ...) {
#'     return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...)/integrate(
#'       Pres_time_adpp, lower = e, upper = s, e = e, s = s,
#'       sp = sp, ...)$value)
#'   }
#'   
#'   lines(x = tt, 
#'         y = Pres_time_adppNorm(tt, s = sim$TS[sp], e = sim$TE[sp], sp = sp),
#'         col="red")
#' }
#' # we provide curves for comparison here, but remember the low sample sizes
#' # (and bins) may affect the quality of the fit
#' 
#' # we can also have a mix of age-independent and age-dependent
#' # models in the same simulation
#' 
#' ###
#' # let us have age-independent sampling before 5my and
#' # age-dependent afterwards
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(2 + 0.1*t)
#' }
#' # note one can also vary the model used for sampling rate
#' # see ?sample.time for examples
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' }
#' 
#' # fixing the extinction times
#' sim$TE[sim$EXTANT] <- 0
#' 
#' # define uniform as above
#' 
#' # preservation function in respect to age
#' # occurrences are uniformly distributed
#' custom.uniform <- function(t, s, e, sp) {
#'   
#'   # make sure it is a valid uniform
#'   if (e >= s) {
#'     message("There is no uniform function with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   res <- dunif(x = t, min = e, max = s)
#'   
#'   return(res)
#' }
#' 
#' # same for TRI
#' 
#' # get the par and par1 vectors
#' # a random quantity
#' par <- runif(n = length(sim$TE), min = sim$TE, max = sim$TS)
#' # its complement to the middle of the lineage's age.
#' # Note that the interaction between these two parameters creates a
#' # deterministic parameter, but inside the function one of them ("par")
#' # is a random parameter
#' par1 <- (((sim$TS - sim$TE)/2) + sim$TE) - par
#' 
#' # preservation function in respect to age, with the "mode" of the triangle
#' # being dependent on par and par1
#' dTRImod <- function(t, s, e, sp) {
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # here is the difference from the function in example 3 and 4
#'   md <- par[sp] + par1[sp]
#'   
#'   # check that md is valid
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t > md & t <= s)
#'   id4 <- which(!(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   res <- vector()
#'   
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e)*(md - e))
#'   res[id2] <- 2 / (s - e)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e)*(s - md))
#'   res[id4] <- 0
#'   
#'   return(res)
#' }
#' 
#' # actual age-dependency defined by a mix
#' dTriAndUniform <- function(t, s, e, sp) {
#'   return(
#'     ifelse(t > 5, custom.uniform(t, s, e, sp),
#'            dTRImod(t, s, e, sp))
#'   )
#' }
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' dt <- sample.clade(sim, rho = rho, tMax = 10, bins = bins,
#'                    adFun = dTriAndUniform, returnTrue = FALSE)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) / 2 + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("t", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[[i]]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)", probability = TRUE)
#'   
#'   # expected curve
#'   tt <- seq(from = sim$TE[sp], to = sim$TS[sp], by = 0.01)
#'   lines(x = tt, y = dTRImod(tt, s = sim$TS[sp], e = sim$TE[sp], sp = sp))
#'   
#'   # need tMax to invert rho
#'   
#'   # getting expected values from age + time dependences:
#'   Pres_time_adpp <- function(t, s, e, sp, ...) {
#'     # correction of scale
#'     rhoMod <- function(t) {
#'       return(rho(tMax - t))
#'     }
#'     return(rhoMod(t)*dTRImod(t = t, s = s, e = e, sp = sp, ...))
#'   }
#'   
#'   # normalizing
#'   Pres_time_adppNorm <- function(t, s, e, sp, ...) {
#'     return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...)/integrate(
#'       Pres_time_adpp, lower = e, upper = s, e = e, s = s,
#'       sp = sp, ...)$value)
#'   }
#'   
#'   lines(x = tt, 
#'         y = Pres_time_adppNorm(tt, s = sim$TS[sp], e = sim$TE[sp], sp = sp),
#'         col="red")
#' }
#' # we provide curves for comparison here, but remember the low sample sizes
#' # (and bins) may affect the quality of the fit
#' 
#' @name sample.clade
#' @rdname sample.clade
#' @export

sample.clade <- function(sim, rho, tMax, S = NULL, envR = NULL, rShifts = NULL,
                         returnTrue = TRUE, returnAll = FALSE, bins = NULL, 
                         adFun = NULL, ...) {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # make S all species if it is NULL
  if (is.null(S)) {
    S <- 1:length(sim$TE)
  }

  # set a default bin
  if (is.null(bins)) {
    bins <- seq(tMax, 0, -0.1)
  } else if (max(bins) < tMax) {
    stop("Bins must include maximum time of simulation")
  }
  
  # if rho is not a constant, apply make.rate to it
  if (!is.numeric(rho) || length(rho) > 1) {
    rho <- make.rate(rho, tMax, envR, rShifts)
  }
  
  # make TE
  sim$TE[sim$EXTANT] <- 0

  # adjusting bins
  bins <- sort(bins, decreasing = TRUE)

  # sample using Poisson process
  
  # independent of age (i.e. occurrences uniformly distributed through the 
  # lineage's age)
  if (is.null(adFun)) { 
    # find occurrences times
    pointEstimates <- sample.time(sim, rho, tMax, S)
    
    # which species left no occurrences
    zeroOccs <- which(lapply(pointEstimates, length) == 0)
    
    # tell the user
    message(paste0(length(zeroOccs), " species left no fossil"))
  } 
  
  # dependent of age (i.e. occurrences distributed through the lineage's age 
  # according to adFun)
  else { 
    # find occurrence times
    pointEstimates <- sample.age(sim = sim, rho = rho, tMax = tMax, S = S,
                                 adFun = adFun, ...)
  }

  # create data frame containing both ranges and time points
  res <- data.frame(matrix(nrow = 0, ncol = 5))
  
  # name the columns
  colnames(res) <- c("Species", "Extant", "SampT", "MaxT", "MinT")
  
  # for each occurrence
  for (i in 1:length(pointEstimates)) {
    
    if (length(pointEstimates[[i]]) > 0) {
      # bin it
      binned_occs <- binner(pointEstimates[[i]], bins = bins)
      
      # counter for point estimates
      counter <- 1
      
      # for each bin
      for (k in 1:(length(bins) - 1)) {
        # if there are occurrences in that bin
        if (binned_occs[k] > 0) {
          # make a row of the data frame
          aux <- data.frame(Species = rep(i, times = binned_occs[k]),
                            Extant = NA, 
                            SampT = pointEstimates[[i]][counter:(counter + binned_occs[k] - 1)],
                            MinT = rep(bins[k + 1], times = binned_occs[k]),
                            MaxT = rep(bins[k], times = binned_occs[k]))
          
          # add row to data frame
          res <- rbind(res, aux)
          
          # increase counter
          counter <- counter + binned_occs[k]
        }
      }
    }
  }
  
  if (nrow(res) > 0) {
    # make the extant column
    res$Extant <- FALSE
    
    # based on the vector in sim
    res$Extant[res$Species %in% which(sim$EXTANT)] <- TRUE
    
    # and the species column
    res$Species <- paste0("t", res$Species)
  }

  if (returnAll) {
    # if returnAll is true, return res as is
    return(res)
  } else if (returnTrue) {
    # if returnTrue is true and returnAll is not, delete age ranges
    return(res[, c("Species", "Extant", "SampT")])
  } else {
    # if neither are true, delete true time
    return(res[, c("Species", "Extant", "MinT", "MaxT")])
  }
}
