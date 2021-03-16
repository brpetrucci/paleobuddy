#' Trait-dependent Birth-Death simulation
#'
#' Simulates a species birth-death process with time and trait dependent rates
#' for any number of starting species. Allows for the speciation/extinction rate
#' to be (1) a constant, or (2) a function of time and a list of trait values.
#' Traits are simulated to evolve under Brownian Motion, Ornstein-Uhlenbeck,
#' Early Burst or MuSSE models (see references). 
#' Allows for constraining results on the number of species at the end of the 
#' simulation, either total or extant. The function can also take an optional 
#' shape argument to generate age-dependence on speciation and/or extinction, 
#' assuming a Weibull distribution as a model of age-dependence. Returns a 
#' \code{sim} object (see \code{?sim}). It may return true extinction times or 
#' simply information on whether species lived after the maximum simulation time,
#' depending on input. For constant rate simulations, see \code{bd.sim.constant}. 
#' For simulations not involving trait dependence, see \code{bd.sim.general}. For
#' a function that unites all scenarios, see \code{bd.sim}. \code{bd.sim} also 
#' allows for extra inputs. For similar flexibility, use \code{make.rate} to 
#' generate the desired rate and then modify it to include trait dependence.
#' Please note while time runs from \code{0} to \code{tMax} in the simulation, it
#' returns speciation/extinction times as \code{tMax} (origin of the group) to 
#' \code{0} (the "present" and end of simulation), so as to conform to other
#' packages in the literature.
#'
#' @inheritParams bd.sim
#' 
#' @param lambda Function to hold the speciation rate over time. It will either be
#' interpreted as an exponential rate, or a Weibull scale if 
#' \code{lShape != NULL}. Can be constant, to allow for mixing of constant and
#' non-constant rates. Can also be a function of a list of traits, 
#' \code{lambda(t, traits)}. For each species a trait evolution simulation will be
#' run for each trait, and then used to calculate the final speciation rate.
#' Note that \code{lambda} should always be greater than or equal to zero.
#'
#' @param mu Similar to above, but for the extinction rate.
#' 
#' Note: rates should be considered as running from \code{0} to \code{tMax}, as
#' the simulation runs in that direction even though the function inverts times
#' before returning in the end.
#' 
#' Note: this function is meant to be called by \code{bd.sim}, so it neither
#' allows for as much flexibility, nor calls \code{make.rate}. If the user wishes
#' to use \code{bd.sim.traits} with environmental or step-function rates, they
#' can generate the rate with \code{make.rate} and supply it to the function.
#' 
#' @param nTraits The number of traits to be considered. \code{lambda} and 
#' \code{mu} need not reference every trait simulated.
#' 
#' @param traitModel The model under which each trait is to be simulated. It can
#' be "BM" for Brownian Motion, "OU" for the Ornstein-Uhlenbeck process, "EB" for
#' the Early Burst process, or "ST" for discrete state trait evolution (MuSSE). 
#' If all traits are supposed to follow the same model, can be just one string.
#' 
#' @param pars A named list of trait parameters. Should be organized as follows. 
#' 
#' \itemize{
#' 
#' \item \code{pars[["BM"]]} should be a named list containing "sigma2" and "X0"
#' fields. See \code{?traits.bm} for more information.
#' 
#' \item \code{pars[["OU"]]} should be a named list containing "sigma2", "X0",
#' "mean" and "theta" fields. See \code{?traits.ou} for more information.
#' 
#' \item \code{pars[["EB"]]} should be a named list containing "sigma2", "X0" and
#' "b" fields. See \code{?traits.eb} for more information.
#' 
#' \item \code{pars[["ST"]]} should be a named list containing "states", "X0" and
#' "Q" fields. Since \code{X0} and \code{Q} are always of length greater than 1,
#' they must be passed inside a list to allow for tests to run correctly. See 
#' \code{?traits.states} for more information.
#' }
#'
#' @return A \code{sim} object, containing extinction times, speciation times,
#' parent, and status information for each species in the simulation. See 
#' \code{?sim}, and a list object with the trait value functions for each
#' species in the simulation.
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @references 
#' 
#' Felsenstein, J. (1973). Maximum-Likelihood estimation of evolutionary trees
#' from continuous chracters. American Journal of Human Genetics. 25(5):471-492.
#'
#' Uhlenbeck, G. E., Ornstein, L. S. (1930). On the theory of Brownian Motion. 
#' Phys. Rev. 36 (5): 823–841.
#' 
#' Felsenstein, J. (1988). Phylogenies and Quantitative Characters. Annual Review
#' of Ecology and Systematics. 19:445-71.
#' 
#' Harmon, L. J., Losos, J. B., Davies, T. D., Gillespie, R. G., Gittleman, J. L.,
#' Jennings, W. B., Hozak, K. H., McPeek, M. A., Moreno-Roark, F., Near, T. J.,
#' Purvis, A., Ricklefs, R. E., Schluter, D., Schulte II, J. A., Seehausen, O.,
#' Sidlauskas, B. L., Torres-Carvajal, O., Weir, J. T., Mooers, A. Ø. (2010).
#' Early Bursts of Body Size and Shape Evolution are Rare in Comparative Data.
#' Evolution 64-8: 2385-2396.
#' 
#' Maddison W.P., Midford P.E., Otto S.P. 2007. Estimating a binary character’s 
#' effect on speciation and extinction. Systematic Biology. 56(5):701.
#' 
#' FitzJohn R.G. 2012. Diversitree: Comparative Phylogenetic Analyses of 
#' Diversification in R. Methods in Ecology and Evolution. 3:1084–1092.
#'
#' @examples
#'
#' ###
#' 
#' @name bd.sim.traits
#' @rdname bd.sim.traits
#' @export

bd.sim.traits <- function(n0, lambda, mu, tMax, 
                          lShape = NULL, mShape = NULL,
                          nTraits = 1, traitModel = "BM",
                          pars = list("BM" = list("sigma2" = 1, "X0" = 0)),
                          nFinal = c(0, Inf), nExtant = c(0, Inf),
                          trueExt = FALSE) {
  # check that n0 is not negative
  if (n0 <= 0) {
    stop("initial number of species must be positive")
  }
  
  # check nFinal's length
  if (length(nFinal) != 2) {
    stop("nFinal must be a vector with a minimum and maximum number 
         of species")
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
  
  # check that traitModel has length 1 or nTraits
  if (length(traitModel) != 1) {
    if (length(traitModel) != nTraits) {
      stop("traitModel must have length 1 or nTraits")
    } 
    ### at some point have to write an error check for pars length
  } else if (nTraits > 1) {
    traitModel <- rep(traitModel, nTraits)
  }
  
  # if it is 1 and nTraits is higher than 1, make it a vector
  
  # initialize test making sure while loop runs
  inBounds <- FALSE
  
  # counter to make sure the nFinal is achievable
  counter <- 1
  
  # make a back up of the rates
  lambda0 <- lambda
  mu0 <- mu
  
  while (!inBounds) {
    # create vectors to hold times of speciation, extinction, 
    # parents and status
    TS <- rep(0, n0)
    TE <- rep(NA, n0)
    parent <- rep(NA, n0)
    isExtant <- rep(TRUE, n0)
    
    # initialize species count
    sCount <- 1

    # initialize list of traits
    traits <- list()
    
    # while we have species to be analyzed still
    while (length(TE) >= sCount) {
      # start at the time of speciation of sCount
      tNow <- TS[sCount]
      
      # if sCount came from another species, it starts
      # its trait evolution with the trait value of that
      # species at time of speciation
      if (!is.na(parent[sCount])) {
        spPars <- pars
        for (m in traitModel) {
          # trait values at given times
          vals <- c()
          
          for (i in 1:nTraits) {
            if (traitModel[i] == m) {
              vals <- c(vals, 
                        traits[[parent[sCount]]][[i]](TS[sCount]))
            }
          }
          
          # assign the corresponding X0 to those values
          spPars[[m]][["X0"]] <- vals
        }
      } else {
        # if it started this simulation, nothing changes
        spPars <- pars
      }

      # run trait evolution and append it to list
      traits[[paste0("t", sCount)]] <- 
        traits.species(tMax = tMax, 
                       tStart = TS[sCount], nTraits = nTraits, 
                       traitModel = traitModel, pars = spPars)
      
      # incorporate trait values in lambda and mu if needed
      
      # if lambda is a function and has more than one argument
      if (!is.numeric(lambda0)) {
        if (length(formals(lambda0)) == 2) {
          # define lambda with trait functions and lambda0
          lambda <- function(t) {
            # for t, get lambda of t and the corresponding trait values
            ifelse(t < TS[sCount], 0,
                   lambda0(t, 
                           unname(unlist(
                             lapply(traits[[sCount]], function(x) x(t))))))
            # need unname and unlist due to how lapply returns lists
            
            # note that if t < TS[sCount] the function has value
            # 0, since the species did not exist
          }
        }

        # make sure it is not negative
        if (optimize(lambda, interval = c(0, 1e10))$objective < 0) {
          stop("speciation rate can never be negative - ensure your function
             is not negative for any trait value")
        }
      }
      
      # if mu is a function and has more than one argument
      if (!is.numeric(mu0)) {
        if (length(formals(mu0)) == 2) {
          # define mu with trait functions and mu0
          mu <- function(t) {
            # for t, get mu of t and the corresponding trait values
            ifelse(t < TS[sCount], 0,
                   mu0(t, 
                           unname(unlist(
                             lapply(traits[[sCount]], function(x) x(t))))))
            # need unname and unlist due to how lapply returns lists
            
            # note that if t < TS[sCount] the function has value
            # 0, since the species did not exist
          }
        }
        
        # make sure it is not negative
        if (optimize(mu, interval = c(0, 1e10))$objective < 0) {
          stop("extinction rate can never be negative - ensure your function
             is not negative for any trait value")
        }
      }

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
  
  res <- list("TRAITS" = traits, "SIM" = sim)
  
  return(res)
}

###
# function to run trait evolution for a given species
traits.species <- function(tMax, tStart, nTraits, traitModel, pars) {
  # create a return list
  traitFuncs <- list()
  
  # modify pars vectors as needed
  
  # make pars vectors the corresponding length if they are constants
  pars <- lapply(c("BM", "OU", "EB", "ST"), function(x) {
    modPars <- pars[[x]]
    for (i in 1:length(modPars)) {
      if (length(modPars[[i]]) != sum(traitModel == x)) {
        if (length(modPars[[i]]) != 1) {
          stop("parameter vectors must be of length 1 or length equal
               to the number of traits under the corresponding model")
        } else {
          modPars[[i]] <- rep(modPars[[i]], sum(traitModel == x))
        }
      }
    }
    return(modPars)
  })
  
  # redo names
  names(pars) <- c("BM", "OU", "EB", "ST")
  
  # take each of them
  bmPars <- pars[["BM"]]
  ouPars <- pars[["OU"]]
  ebPars <- pars[["EB"]]
  stPars <- pars[["ST"]]
  
  # start a count to see where we are on the corresponding model counts
  bmCount <- 1
  ouCount <- 1
  ebCount <- 1
  stCount <- 1
  
  for (i in 1:nTraits) {
    if (traitModel[i] == "BM") {
      # take the (bmCount)th pars and run evolution
      traitFuncs[paste0("trait", i)] <- 
        traits.bm(tMax = tMax, nTraits = 1, tStart = tStart, 
                  X0 = bmPars[["X0"]][bmCount],
                  sigma2 = bmPars[["sigma2"]][bmCount],
                  bounds = bmPars[["bounds"]][[bmCount]])
      
      # increase the count
      bmCount <- bmCount + 1
    }
    
    else if (traitModel[i] == "OU") {
      # take the (ouCount)th pars and run evolution
      traitFuncs[paste0("trait", i)] <-
        traits.ou(tMax = tMax, nTraits = 1, tStart = tStart,
                  sigma2 = ouPars[["sigma2"]][ouCount], 
                  theta = ouPars[["theta"]][ouCount], 
                  mean = ouPars[["mean"]][ouCount], X0 = ouPars[["X0"]][ouCount])
      
      # increase the count
      ouCount <- ouCount + 1
    }
    
    else if (traitModel[i] == "EB") {
      # take the (ebCount)th pars and run evolution
      traitFuncs[paste0("trait", i)] <-
        traits.eb(tMax = tMax, tStart = tStart, nTraits = 1, 
                  sigma2 = ebPars[["sigma2"]][ebCount], 
                  b = ebPars[["b"]][ebCount], X0 = ebPars[["X0"]][ebCount])
      
      # increase the count
      ebCount <- ebCount + 1
    }
    
    else if (traitModel[i] == "ST") {
      # take the (stCount)th pars and run evolution
      traitFuncs[paste0("trait", i)] <-
        traits.states(tMax = tMax, tStart = tStart, nTraits = 1, 
                     states = stPars[["states"]][[stCount]],
                     Q = stPars[["Q"]][[stCount]], 
                     X0 = stPars[["X0"]][stCount])
      
      # increase the count
      stCount <- stCount + 1
    }
  }
  
  return(traitFuncs)
}