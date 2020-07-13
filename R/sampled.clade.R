#' General rate species sampling
#'
#' \code{sample.clade} takes times of speciation and extinction, information to
#' create a sampling rate with \code{make.rate}, a vector of geologic time 
#' intervals, and whether one wants the true return times or a range based on 
#' \code{bins}, and returns a data frame with minimum and maximum times for
#' occurrences (\code{returnTrue == FALSE}) or exact occurrence times 
#' (\code{returnTrue == TRUE}) for each species.
#'
#' @param S a list species numbers to be sampled. Could be only a subset of the
#' species if the user wishes. The default is all species in \code{sim}.
#'
#' @param sim a simulation, usually the output of \code{bd.sim}.
#'
#' @param bins a vector of time intervals corresponding to geological time
#' ranges. If \code{returnTrue} is false, \code{sample.clade} returns the
#' occurrence times as ranges. In this way, we simulate the granularity in
#' real world fossil records. If \code{returnTrue} is true, this is ignored.
#'
#' @param rr a sampling rate function. May be a constant, a time-dependent
#' function, a function dependent on time and environment, or a vector of
#' rates corresponding to the times in \code{rShifts}.
#' 
#' Note: must be a constant if \code{dFun} is not \code{NULL}.
#'
#' @param tMax the maximum simulation time, used by \code{rexp.var}.
#' @param envRR a matrix containing time points and values of an enviromental
#' variable, like temperature, for each time point. This will be used to create
#' a sampling rate, so \code{rr} must be a function of time and said variable
#' if \code{envRR} is not NULL.
#'
#' @param rShifts vector of rate shifts. First element must be the sstarting
#' time for the simulation (0 or tMax). It must have the same length as
#' \code{rr}. E.g. \code{rr = c(0.1, 0.2, 0.1)}, \code{rShifts = c(0, 10, 20)}
#' means the sampling rate will be 0.1 from 0 to 10, 0.2 from 10 to 20, and 0.1
#' from 20 to \code{tMax}. It would also be identical, in this case, to use
#' \code{pshifts = c(tMax, tMax - 10, tMax - 20)}.
#'
#' Note that using this method for step-function rates is currently slower than
#' using \code{ifelse}.
#'
#' @param returnTrue if set to \code{TRUE}, the returned data frame will
#' contain true times of sampling. If set to \code{FALSE}, we call 
#' \code{binner} and the returned data frame will contain ranges of sampling
#' times based on \code{bins}.
#'
#' @param dFun A density function representing the age-dependent
#' preservation model. It must be a density function, and consequently
#'
#' \describe{
#'
#' \item{1.}{integrate to 1 (though this condition is not verified by the 
#' function, it is the user's responsibility to check this property)}
#'
#' \item{2.}{describe the density of sampling a lineage in a given point \code{t}
#' in geological time}
#'
#' \item{3.}{be parametrized in absolute geological time (i.e. should be relative
#' to absolute geological time, in Mya)}
#'
#' \item{4.}{should be limited between \code{s} (i.e. the lineage's 
#' speciation/origination geological time) and \code{e} (i.e. the lineage's 
#' extinction geological time), with \code{s} > \code{e}}
#'
#' \item{5.}{include the arguments \code{t}, \code{s}, \code{e} and \code{sp}}}
#'
#' @param dFunMax a function that calculates the maximum (density) value
#' of \code{dFun} using its arguments. It can also be a number representing the
#' maximum density.
#'
#' Note that if it is not provided, it will be approximated numerically, leading 
#' to longer running times.
#'
#' @param ... additional parameters related to \code{dFun} and \code{dFunMax}.
#'
#' @return a data frames containing species names/numbers, whether each species
#' is extant, and either the true occurrence times of species or a range of
#' occurrence times based on \code{bins}.
#'
#' @author written by Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @examples
#'
#' # in all the coming examples we show histograms for visualization, but note
#' # they might not look exactly as intended given the low sample size. See
#' # vignettes for more thorough testing
#' 
#' ###
#' # we can start with a constant case
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # sampling rate
#' r <- 2
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' # find the occurrence data frame
#' dt <- sample.clade(1:length(sim$TE), sim, r, tMax = 10, bins = bins)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("spp_", "", ids[i])))
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
#' # sampling can be any function of time in the non-age dependent case, of course
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # sampling rate
#' r <- function(t) {
#'   return(3 - 0.25*t)
#' }
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' # find the occurrence data frame
#' dt <- sample.clade(1:length(sim$TE), sim, r, tMax = 10, bins = bins)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("spp_", "", ids[i])))
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
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we will use the less efficient method of creating a step function
#' # one could instead use ifelse()
#' 
#' # rates vector
#' rlist <- c(1, 2, 0.5)
#' 
#' # rate shifts vector
#' rShifts <- c(0, 4, 8)
#' 
#' # make it a function so we can plot it
#' r <- make.rate(rlist, 10, fShifts=rShifts)
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' # find the occurrence data frame
#' dt <- sample.clade(1:length(sim$TE), sim, r, tMax = 10, bins = bins)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("spp_", "", ids[i])))
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
#' # finally, \code{sample.clade} also accepts an environmental variable
#' if (requireNamespace("RPANDA", quietly = TRUE)) {
#'   # get temperature data
#'   data(InfTemp, package = "RPANDA")
#'   
#'   # simulate a group
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#'   
#'   # in case first simulation is short-lived
#'   while ((sim$TS[1] - sim$TE[1]) < 10) {
#'     sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#'   }
#'   
#'   # make temperature the environmental dependency of r
#'   envR <- InfTemp
#'   
#'   # we can then make sampling dependent on the temperature
#'   r <- function(t, env) {
#'     return(0.5*env)
#'   }
#'   
#'   # make it a function so we can plot it
#'   rr <- make.rate(r, envF = envR)
#'   
#'   # let us check that r is high enough to see a pattern
#'   plot(1:10, rr(1:10), type = 'l', main = "Sampling rate",
#'        xlab = "My", ylab = "r")
#'   
#'   # the resolution of the fossil dataset:
#'   bins <- seq(from = 10, to = 0,
#'               by = -0.1)
#'   # note that we will provide a very high resolution to test the function
#'   
#'   # find the occurrence data frame
#'   dt <- sample.clade(1:length(sim$TE), sim, r, tMax = 10, envRR = envR,
#'                     bins = bins)
#'   
#'   # extract species identity
#'   ids <- unique(dt$Species)
#'   
#'   # approximate sampling time (since it is a range)
#'   mids <- (dt$MaxT - dt$MinT) + dt$MinT
#'   
#'   # for each species
#'   for (i in 1:length(ids)) {
#'     # get the species number
#'     sp <- unique(as.numeric(gsub("spp_", "", ids[i])))
#'     
#'     # check the histogram
#'     hist(mids[dt$Species == ids[i]],
#'          main = paste0("spp = ", sp, "; duration ~ ",
#'                        round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'          xlab = "Time (My)",
#'          xlim = c(sim$TS[i], sim$TE[i]))
#'   }
#' }
#' 
#' # we will now do some tests with age-dependent rates. For more details,
#' # check sample.adpp.
#' 
#' ###
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
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
#' # function to calculate max of the PERT
#' dPERTmax <- function(s, e, sp) {
#'   return(((s - e) / 2) + e)
#' }
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' dt <- sample.clade(1:length(sim$TE), sim, rr = 3, tMax = 10, bins = bins,
#'                   dFun = dPERT, dFunMax = dPERTmax)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("spp_", "", ids[i])))
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
#' # (and bins) may affect the quality of the fit. See vignettes for more 
#' # thorough testing
#' 
#' ###
#' # now, a hat-shaped increase through the duration of a species dependent on two
#' # parameters
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - sim$TE[1]) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # preservation function in respect to age, with the "mode" of the triangle
#' # being exactly at the last quarter of the duration of EACH lineage.
#' dTRImod2<-function(t, s, e, sp) {
#'   
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
#'   #for more details in this function, see example 3 and 4
#' }
#' 
#' # maximum function
#' dTRImaxmod2 <- function(s, e, sp) {
#'   return(2 / (s - e))
#' }
#' 
#' # a random point inside each lineage's duration
#' par <- runif (n = length(sim$TE), min = sim$TE, max = sim$TS)
#' 
#' # a distance between "par" and the lineage's duration middle
#' par1 <- (((sim$TS - sim$TE) / 2) + sim$TE) - par
#' 
#' # the resolution of the fossil dataset:
#' bins <- seq(from = 10, to = 0,
#'             by = -0.1)
#' # note that we will provide a very high resolution to test the function
#' 
#' dt <- sample.clade(1:length(sim$TE), sim, rr = 4, tMax = 10, bins = bins,
#'                   dFun = dTRImod2, dFunMax = dTRImaxmod2)
#' 
#' # extract species identity
#' ids <- unique(dt$Species)
#' 
#' # approximate sampling time (since it is a range)
#' mids <- (dt$MaxT - dt$MinT) + dt$MinT
#' 
#' # for each species
#' for (i in 1:length(ids)) {
#'   # get the species number
#'   sp <- unique(as.numeric(gsub("spp_", "", ids[i])))
#'   
#'   # check the histogram
#'   hist(mids[dt$Species == ids[[i]]],
#'        main = paste0("spp = ", sp, "; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
#'        xlab = "Time (My)", probability = TRUE)
#'   
#'   # expected curve
#'   curve(dTRImod2(x, e=sim$TE[sp], s=sim$TS[sp], sp=sp),from = sim$TE[sp],
#'         to = sim$TS[sp], add=TRUE, col="red", n = 100)
#' }
#' # we provide curves for comparison here, but remember the low sample sizes 
#' # (and bins) may affect the quality of the fit. See vignettes for more 
#' # thorough testing
#' 
#' @name sample.clade
#' @rdname sample.clade
#' @export

sample.clade <- function(S = NULL, sim, rr, tMax, envRR = NULL, rShifts = NULL,
                         returnTrue = FALSE, bins = NULL, 
                         dFun = NULL, dFunMax = NULL,...) {
  # make S all species if it is NULL
  if (is.null(S)) {
    S = 1:length(sim$TE)
  }
  
  # get the speciation and extinction times vectors
  TE <- sim$TE[S]
  TS <- sim$TS[S]

  # check if it is age-dependent
  if (is.null(dFun)) {
    rr <- make.rate(rr, tMax, envRR, rShifts)
  } else {
    if (!is.numeric(rr) | length(rr)>1) {
      stop("ADPP cannot be used with time-varing preservation rates")
    }
  }
  
  # check if we have bins if we need them
  if (!returnTrue & is.null(bins)) {
    stop("sample.clade needs a bins vector to returned binned samples")
  }

  # adjusting bins
  bins <- sort(bins, decreasing = TRUE)

  # sample using Poisson process
  
  # independent of age (i.e. occurrences uniformly distributed through the 
  # lineage's age)
  if (is.null(dFun)) { 
    # find occurrences times
    pointEstimates <- lapply(S, sample, TE = TE, TS = TS, rr = rr, tMax = tMax)
    
    # which species left no occurrences
    zeroOccs <- which(lapply(pointEstimates, length) == 0)
    
    # tell the user
    message(paste0(length(zeroOccs), " species left no fossil"))
  } 
  
  #dependent of age (i.e. occurrences distributed through the lineage's age 
  # according to dFun)
  else { 
    # find occurrence times
    pointEstimates <- sample.adpp(S, TS = TS, TE = TE, rr = rr, 
                                 dFun = dFun, dFunMax = dFunMax, ...)
  }

  # wrapping data

  # output as fossil occurrence binned within bins/bins
  if (!returnTrue) {
    # create the data frame
    res <- data.frame(matrix(nrow = 0, ncol = 4))
    
    # name the columns
    colnames(res) <- c("Species", "Extant", "MaxT", "MinT")
    
    # for each occurrence
    for (i in 1:length(pointEstimates)) {
      
      # bin it
      binned_occs <- binner(pointEstimates[[i]], bins = bins)
      
      # for each bin
      for (k in 1:(length(bins) - 1)) {
        # if there are occurrences in that bin
        if (binned_occs[k] > 0) {
          # make a row of the data frame
          aux <- data.frame(Species = i, 
                            Extant = NA, 
                            MaxT = rep(bins[k], times = binned_occs[k]), 
                            MinT = bins[k + 1])
          
          # add row to data frame
          res <- rbind(res, aux)
        }
      }
    }
    # make the extant column
    res$Extant <- FALSE
    
    # based on the list in sim
    res$Extant[res$Species %in% which(sim$EXTANT)] <- TRUE

    # and the species column
    res$Species <- paste0("spp_", res$Species)
  } 
  
  # if user wants the true times of occurrences, get a a data frame with the 
  # real sampling times only
  else {
    # create the data frame - one less column
    res <- data.frame(matrix(nrow = length(unlist(pointEstimates)), ncol = 3))
    
    # name the columns
    colnames(res) <- c("Species", "Extant", "SampT")
    
    # if there are occurrences
    if (nrow(res) > 1) {
      # make species column
      res$Species <- rep(S, times = lapply(pointEstimates, length))
      
      # make the extant column
      res$Extant <- FALSE
      
      # based on the list in sim
      res$Extant[res$Species %in% which(sim$EXTANT)] <- TRUE
      
      # name the species 
      res$Species <- paste0("spp_", res$Species)
      
      # make the sampling times column
      res$SampT <- unlist(pointEstimates)
    }
  }

  return(res)
}
