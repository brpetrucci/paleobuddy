#' Time-dependent and age-dependent rate species sampling
#' 
#' Generates a vector of occurrence times for a species in a simulation using 
#' a Poisson process coupled with age-dependent fossil sampling. Allows for the 
#' Poisson rate to be (1) constant or (2) time-dependent, and additionally 
#' allows the distribution of occurrences throughout a species age to be any 
#' function between speciation and extinction time. For more flexibility, see 
#' \code{make.rate} and \code{sample.clade}. Note that while the simulation 
#' occurs in forward time, we return (both in birth-death functions and here)
#' results in backwards time, so that time is inverted using \code{tMax} both 
#' at the beginning and end of the function.
#' 
#' @param rho Average fossil sampling rate (per species per million years) over 
#' time. It can be a \code{numeric} describing a constant rate or a 
#' \code{function(t)} describing the variation in sampling over time. For more 
#' flexibility on sampling, see \code{make.rate} to create more complex rates.
#' If \code{adFun} is supplied, it will be used to find the number of 
#' occurrences during the species duration, and \code{adFun} will determine 
#' their distribution. Note that \code{rho} should always be greater than or 
#' equal to zero.
#'
#' @param adFun A density function representing the age-dependent preservation
#' model. It must be a density function, and consequently integrate to 1 
#' (though this condition is not verified by the function). If not provided, a 
#' uniform distribution will be used by default. The function must also 
#' have the following properties:
#'
#' \itemize{
#'
#' \item Return a vector of preservation densities for each time in a given
#' vector \code{t} in geological time. 
#'
#' \item Be parametrized in absolute geological time (i.e. should be relative 
#' to absolute geological time, in Mya, \emph{not} the lineage's age). Because 
#' of this, it is assumed to go from \code{tMax} to \code{0}, as opposed to 
#' most functions in the package.
#'
#' \item Should be limited between \code{s} (i.e. the lineage's 
#' speciation/birth) and \code{e} (i.e. the lineage's extinction/death), 
#' with \code{s} > \code{e}.
#'
#' \item Include the arguments \code{t}, \code{s}, \code{e} and \code{sp}. 
#' The argument sp is used to pass species-specific parameters (see examples),
#' allowing for \code{dFun} to be species-inhomogeneous.
#' }
#'
#' @param ... Additional parameters used by the user's function (i.e.,
#' \code{adFun})
#' 
#' @inheritParams sample.time
#'
#' @return A list of occurrences for that species, given the time, age 
#' and species-specific conditions assigned by the user.
#'
#' @author Matheus Januario
#'
#' @examples
#' 
#Please note the shape of the fossilization becomes easier to see when rho is very high (e.g., >100). We will not do this here due to CRAN constraints on example length, but users are encoraged to try. Histograms can be used to visualize the sampling though age, as the examples do in the 2nd example (PERT function) - we will not do this is all examples to get a smaller number of code lines as examples
#' 
#' # vector of times
#' time <- seq(10, 0, -0.1)
#' 
#' ###
#' # we can start with a constant Poisson rate
#' 
#' # sampling rate
#' rho <- 3
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # add a hat-shaped increase through the duration of a species
#' 
#' # here we will use the PERT function as described
#' # in Silvestro et al 2014
#' 
#' # age-dependence distribution
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
#' # plot how this functions translate preservation and age for an example species:
#' S_sp <- 10 #speciation time
#' E_sp <- 5 #extinction time
#' 
#' plot(time, rev(dPERT(t = time, s = S_sp, e = E_sp, a = 1)), main = "Age-dependence distribution",
#'      xlab = "Species age (My)", ylab = "Density", 
#'      xlim = c(0, S_sp-E_sp), type = "l")
#' 
#' # sample first and third species only
#' occs <- sample.age(sim, rho = 3, tMax=10, dPERT, S = c(1,3))
#' 
#' #making histograms of each sp occs:
#' par(mfrow=c(1,2))
#' hist(occs[[1]], main="Occs of sp1")
#' hist(occs[[3]], main="Occs of sp3")
#' par(mfrow=c(1,1))
#' 
#' ###
#' # now we can try time-dependent, age-indenpendent sampling
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(1 + 0.5*t) 
#' }
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # age-independence distribution (a uniform density distribution in age) inputted in the format that the function needs
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
#' # this is the same as not supplying adFun, just an illustration
#' 
#' # sample first 3 species
#' occs <- sample.age(sim, rho, tMax = 10, S = 1:3,
#'                    adFun = custom.uniform)
#' 
#' par(mfrow=c(1,3))
#' hist(occs[[1]], main="Occs of sp1")
#' hist(occs[[2]], main="Occs of sp2")
#' hist(occs[[3]], main="Occs of sp3")
#' par(mfrow=c(1,1))
#' 
#' ###
#' # we can have more parameters on adFun
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(1 + 0.5*t) 
#' }
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # here we will use the triangular distribution. We have some empirical 
#' # evidence that taxa occurrences might present a triangular shape,
#' # see Zliobaite et al 2017
#' 
#' # age-dependence distribution
#' dTRI <- function(t, s, e, sp, md) {
#'   
#'   # please note ths function is inverted. The correspondence would be:
#'   # s = b maximum
#'   # e = a minimum
#'   # md = c distribution's mode
#'   
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # another condition we must check
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # needed to vectorize the function:
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t > md & t <= s)
#'   id4 <- which( !(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   # actually vetorizing function
#'   res <- vector()
#'   
#'   # (t >= e & t < md)
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e) * (md - e))
#'   
#'   #(t == md)
#'   res[id2] <- 2 / (s - e)
#'   
#'   #(md < t & t <= s)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e) * (s - md))
#'   
#'   #outside function's limits
#'   res[id4] <- 0
#'   
#'   return(res)
#' }
#' 
#' # plot for the first species to take a look
#' S_sp <- 10 #speciation time
#' E_sp <- 5 #extinction time
#' age=seq(S_sp, E_sp, by=-0.1)
#' 
#' plot(time, rev(dTRI(time, S_sp, E_sp, 1, 9)), 
#'      main = "Age-dependence distribution",
#'      xlab = "Species age (My)", ylab = "Density", 
#'      xlim = c(0, 10), type = "l")
#' 
#' # sample first species
#' occs <- sample.age(sim = sim, rho = rho, adFun = dTRI, 
#'                    tMax = 10, S = 1, md = 8)
#' # note we are providing the mode (in absolute geologic time, not in "age" - for this ecxample see below) for the triangular sampling as an argument
#' 
#' ###
#' # we can also have a hat-shaped increase through the duration of a species
#' # with more parameters than TS and TE, but with the parameters relating to
#' # the relative age of each lineage
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(1 + 0.5*t) 
#' }
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # age-dependence distribution, with the "mde" of the triangle 
#' # being exactly at the last quarter of the duration of EACH lineage
#' dTRImod1 <- function(t, s, e, sp) {
#'   # note that now we don't have the "md" parameter here,
#'   # but it is calculated inside the function
#'   
#'   # check if it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # calculate md
#'   md <- ((s - e) / 4) + e
#'   # md is at the last quarter of the duration of the lineage
#'   
#'   # please note that the same logic can be used to sample parameters
#'   # internally in the function, running for instance:
#'   # md <- runif (n = 1, min = e, max = s)
#'   
#'   # check it is a valid md
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # needed to vectorize function
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t>md & t <= s)
#'   id4 <- which( !(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   # vectorize the function
#'   res<-vector()
#'   
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e) * (md - e))
#'   res[id2] <- 2 / (s - e)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e) * (s - md))
#'   res[id4] <- 0
#'   
#'   return(res)
#' }
#' 
#' # plot for the first species to take a look
#' plot(time, rev(dTRImod1(time, 10, 0, 1)), 
#'      main = "Age-dependence distribution",
#'      xlab = "Species age (My)", ylab = "Density", 
#'      xlim = c(0, 10), type = "l")
#' 
#' # sample first two species
#' occs <- sample.age(sim = sim, rho = rho, adFun = dTRImod1, tMax = 10,
#'                    S = 1:2)
#' 
#' # please note in this case, function dTRImod1 fixes "md" at the last quarter
#' # of the duration of the lineage
#' 
#' ###
#' # the parameters of adFun can also relate to each specific lineage,
#' # which is useful when using variable parameters for each species
#' # to keep track of those parameters after the sampling is over
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(1 + 0.5*t) 
#' }
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # get the par and par1 vectors
#' 
#' # a random quantity
#' par <- runif(n = length(sim$TE), min = 0, max = sim$TS) 
#' # its complement to the middle of the lineage's age. 
#' par1 <- sim$TS / 2 - par 
#' # note that the interaction between these two parameters creates a 
#' # deterministic parameter, but inside the function one of them ("par") 
#' # is a random parameter
#' 
#' # age-dependence distribution, with the "mode" of the triangle
#' # being dependent on par and par1
#' dTRImod2 <- function(t, s, e, sp) {
#'   
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # md depends on parameters
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
#' # plot for the first species to take a look
#' plot(time, rev(dTRImod2(time, 10, 0, 1)), 
#'      main = "Age-dependence distribution",
#'      xlab = "Species age (My)", ylab = "Density", 
#'      xlim = c(0, 10), type = "l")
#' 
#' # sample the first 3 species
#' occs <- sample.age(sim = sim, rho = rho, adFun = dTRImod2, tMax = 10,
#'                    S = 1:3)
#' 
#' # we can also have a mix of age-independent and age-dependent
#' # models in the same simulation
#' 
#' ###
#' # we can have age-independent sampling before 5my and
#' # age-dependent afterwards
#' 
#' # sampling rate
#' rho <- function(t) {
#'   return(1 + 0.5*t)
#' }
#' # note one can also vary the model used for sampling rate
#' # see ?sample.time for examples
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10)
#' 
#' # define uniform as above
#' 
#' # age-dependent distribution
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
#' 
#' # a random quantity
#' par <- runif(n = length(sim$TE), min = 0, max = sim$TS) 
#' # its complement to the middle of the lineage's age. 
#' par1 <- sim$TS / 2 - par 
#' # note that the interaction between these two parameters creates a 
#' # deterministic parameter, but inside the function one of them ("par") 
#' # is a random parameter
#' 
#' # age-dependence distribution, with the "mode" of the triangle
#' # being dependent on par and par1
#' dTRImod2 <- function(t, s, e, sp) {
#'   
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # md depends on parameters
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
#'            dTRImod2(t, s, e, sp))
#'   )
#' }
#' 
#' # plot for the first species to take a look
#' plot(time, rev(dTriAndUniform(time, 10, 0, 1)), 
#'      main = "Age-dependence distribution",
#'      xlab = "Species age (My)", ylab = "Density", 
#'      xlim = c(0, 10), type = "l")
#' 
#' # sample the first species
#' occs <- sample.age(sim = sim, rho = rho, adFun = dTriAndUniform,
#'                    tMax = 10, S = 1)
#' 
#' # please note in this case, function dTRImod2 fixes "md" at the last quarter
#' # of the duration of the lineage

sample.age.time <- function(sim, rho, tMax, S = NULL, adFun = NULL, ...){
  # checking input
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # convert NA elements of TE to 0
  # those automatically to zeor and print warning
  if (sum(is.na(sim$TE)) > 0) {
    message("TEs of extant species will be converted to 0")
    sim$TE[sim$EXTANT] <- 0
  }
  
  # if rho is constant, make it a function
  if (is.numeric(rho)) {
    r <- rho
    rho <- Vectorize(function(t) r)
  }
  
  # if S is not given
  if (is.null(S)) {
    S <- 1:length(sim$TE)
  }
  
  # get extinction and speciation times
  TE <- sim$TE
  TS <- sim$TS
  
  # check if adFun is a function with the right parameters
  if (!(is.null(adFun))) {
    if (sum(c("t", "s", "e", "sp") %in% names(formals(adFun))) < 3) {
      stop("adFun must have \"t\", \"s\", \"e\", and  \"sp\" parameters")
    }
  }
  
  # if adFun is null, it is the same as sample.time
  if (is.null(adFun)) {
    message("Preservation will be Age-independent \n")
    
    # run sample.time
    occs <- sample.time(sim, rho, tMax, S)
    
    # check how many species left no occurrences
    zeroOccs <- which(lapply(occs, length) == 0)
    message(paste0(length(zeroOccs), " species left no fossil"))
    return(occs)
  }
  
  # create function combining adFun and rho
  Pres_time_adpp <- function(t, s, e, sp, ...) {
    # this correction is needed because "t" in rho means
    # "time since clade origin" and NOT absolute geologic time in Mya
    rhoMod <- function(t){ 
      return(rho(tMax - t))
    }
    return(rhoMod(t)*adFun(t = t, s = s, e = e, sp = sp, ...))
  }
  
  # normalize it
  Pres_time_adppNorm <- function(t, s, e, sp,...) {
    return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...) / 
             integrate(Pres_time_adpp, lower = e, upper = s, s = s, 
                       e = e, sp = sp, ...)$value)
  }
  
  # get the number of occurrences following a time-varying poisson process
  Noccs <- vector(length = length(sim$TE))
  Noccs <- unlist(lapply(sample.time(sim = sim, rho = rho, tMax = tMax,
                                     S = S), FUN = length))
  
  occs <- list()
  
  # for each species to sample
  for (sp in S) {
    # get the maximum of function for that specific species
    adFunMaxApprox <- optimize(Pres_time_adppNorm, interval = c(TE[sp],TS[sp]),
                               e = TE[sp], s = TS[sp], sp = sp, ..., 
                               maximum = T)$maximum
    
    nOccs <- Noccs[sp]
    
    # place occurrences along lineage's age using the accept-reject method
    res <- vector()
    while (length(res) < nOccs) {
      t <- runif(n = 1, max = TS[sp], min = TE[sp])
      
      test <- runif(n = 1, min = 0, max = adFunMaxApprox*1.05) <= 
        Pres_time_adppNorm(t = t, s = TS[sp], e = TE[sp], sp = sp, ...)
      
      test <- runif(n = 1, min = 0, max = Pres_time_adppNorm(t = adFunMaxApprox, 
          s = TS[sp], e = TE[sp], sp = sp, ...)*1.05) <= 
        Pres_time_adppNorm(t = t,
            s = TS[sp], e = TE[sp], sp = sp, ...)
      
      if (test) {
        res <- c(res, t)
      }
    }
    occs[[sp]] <- res
  }
  
  zeroOccs <- which(lapply(occs, length) == 0)
  message(paste0(length(zeroOccs), " species left no fossil"))
  return(occs)
}
