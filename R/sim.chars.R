#' Morphological and molecular character simulation
#' 
#' Description
#' 
#' #at param...
#' 
#' @author Bruno do Rosario Petrucci
#' 
#' #at references...
#' 
#' #at examples...
#' 
#' @name sim.chars
#' @rdname sim.chars
#' @export

sim.chars <- function(sim, sample, nTraits = NULL, Q, 
                      nStates, X0, species_clock,
                      var_min = 0.005, labels = NULL, sim_max = NULL) {
  #### TODO: ERROR CHECKS ####
  # # check that all parameters are the necessary dimensions
  # if (!(length(nStates) %in% c(1, nTraits)) ||
  #     !(length(X0) %in% c(1, nTraits))) {
  #   stop("nStates and X0 must be either of length 1 or nTraits")
  # }
  # if (length(nStates) == 1) nStates <- rep(nStates, nTraits)
  # if (length(X0) == 1) X0 <- rep(X0, nTraits)
  # if (!(species_clock %in% c(1, length(sim$TS)))) {
  #   stop("species_clock must be either of length 1 or the number of species in sim")
  # }
  # if (length(species_clock) == 1) species_clock <- rep(species_clock,
  #                                                      length(sim$TS))
  # if (typeof(Q) == "closure") {
  #   for (i in 1:length(nStates)) {
  #     # get this Q
  #     Q_i <- Q(nStates[i])
  #     
  #     # check that it is the right type and dimensions
  #     if (!(typeof(Q_i) == "double") &&
  #         (nrow(Q_i) == nStates[i] && ncol(Q_i) == nStates[i] && 
  #          all(diag(Q_i) == 0))) {
  #       stop("Q must return a 0-diagonal matrix of dimensions NxN when applied to N")
  #     }
  #   }
  # } else {
  #   if (!(typeof(Q) == "double") && (length(unique(nStates)) == 1) &&
  #       (nrow(Q) == nStates[1] && ncol(Q) == nStates[1] && 
  #        all(diag(Q) == 0))) {
  #     stop("Q must have dimensions NxN, where N is nStates")
  #   }
  # }
  
  # figure out species for which we need to run simulation
  all_species <- paste0("t", 1:length(sim$TS))
  sampled_species <- all_species[sim$EXTANT | all_species %in% sample$Species]
  
  # create chars matrix
  chars <- data.frame(matrix(nrow = length(sampled_species), 
                             ncol = 1 + ifelse(is.null(nTraits),
                                               sum(nStates), nTraits)))
  chars[, 1] <- sampled_species
  
  # check whether we have a fixed number of traits
  if (is.null(nTraits)) {
    # start trait counter
    trait_count <- 1
    
    # if not, iterate through nStates -- assume first field is 0
    # (can change that later to make it easy to get invariant characters)
    for (i in 2:length(nStates)) {
      # how many traits do we want with i states?
      n_traits_i <- nStates[i]
      
      # if none, next
      if (n_traits_i == 0) next
      
      # get Q matrix 
      Q_i <- Q(i)
      
      # start sim counter--we don't want to try to run this for too long
      sim_count <- 1
      
      # iterate through  number of traits
      for (j in 1:n_traits_i) {
        # initialize condition to false
        cond <- FALSE
        
        # until we have the number we want
        while (!cond) {
          # get this character 
          char_n <- sim.char.one(sim, sample, i, X0, Q_i, 
                                 species_clock, labels)
          # allow for invariant character
          
          # increase sim_count
          sim_count <- sim_count + 1
          if (!is.null(sim_max) && sim_count > sim_max) return(NULL)
          
          # condition on minimum variation
          cond <- all(sapply(0:(i - 1), function(x) sum(char_n[, 2] == x)) / 
                        nrow(char_n) > var_min)
        }
        
        # cbind it to chars
        chars[, 1 + trait_count] <- char_n$Char
        trait_count <- trait_count + 1
      }
    }
  } else {
    # if nStates is of length 1, make it the length of nTraits
    # same for X0
    if (length(nStates) == 1) nStates <- rep(nStates, nTraits)
    if (length(X0) == 1) X0 <- rep(X0, nTraits)
    
    # iterate through the number of traits
    for (i in 1:nTraits) {
      # get rate matrix
      Q_i <- Q(nStates[i])
      
      # get this character 
      char_n <- sim.char.one(sim, sample, nStates[i], X0[i], Q_i,
                             species_clock, labels)
      
      # cbind it to chars
      chars[, 1 + i] <- char_n$Char
    }
  }
  
  # make rownames species and delete first column
  rownames(chars) <- chars[, 1]
  chars <- chars[, -1]
  
  # make colnames just trait numbers
  colnames(chars) <- 1:ncol(chars)
  
  # return chars
  return(chars)
}

sim.char.one <- function(sim, sample, nStates, X0, Q, species_clock, 
                         labels = NULL) {
  # check if Q has the correct number of dimensions
  if (nrow(Q) != nStates || nrow(Q) != nStates) {
    # error
    stop("Make sure Q is square and has row/column numbers equal to nStates")
  }
  # check that species_clock is the length of 1 or number of species
  if (length(species_clock) != length(sim$TS)) {
    if (length(species_clock) != 1) {
      stop("species_clock must be of length 1 or the number of species in sim")
    } else {
      species_clock <- rep(species_clock, length(sim$TS))
    }
  }
  
  # get maximum simulation time
  tMax <- max(sim$TS)
  
  # invert sim and sample times since we're simulating traits forward
  sim$TS <- tMax - sim$TS
  sim$TE <- tMax - sim$TE
  sample$SampT <- tMax - sample$SampT
  
  # make a return data frame
  char <- data.frame(matrix(nrow = 0, ncol = 2))
  
  # create a list for the trait values of each species
  traits <- vector("list", length(sim$TS))
  
  # iterate through species in sim
  for (i in 1:length(sim$TS)) {
    # find initial value
    X0_sp <- ifelse(is.na(sim$PAR[i]), X0,
                    traits[[sim$PAR[i]]]$value[
                      findInterval(sim$TS[i], traits[[sim$PAR[i]]]$min)
                    ])
    
    # get time to simulate to
    tMax_sp <- ifelse(is.na(sim$TE[i]), tMax, sim$TE[i])
    
    # get trait values for this species
    traits[[i]] <- sim.char.species(tMax_sp, sim$TS[i], nStates, X0_sp,
                                    Q * species_clock[i])
    
    # get samples from this species
    sample_sp <- sample[sample$Species == paste0("t", i), ]
    
    # if species is extant, get final state
    if (sim$EXTANT[i]) {
      char <- rbind(char, c(paste0("t", i), tail(traits[[i]]$value, 1)))
    } else if (nrow(sample_sp) > 0) {
      # if not, choose a random time from within fossil range
      char_samp_time <- runif(1, min(sample_sp$SampT), max(sample_sp$SampT))
      
      # and sample state at that time
      char <- rbind(char, c(paste0("t", i),
                            traits[[i]]$value[findInterval(char_samp_time,
                                                           traits[[i]]$min)]))
    }
  }

  # name char
  colnames(char) <- c("Species", "Char")

  # if labels is not null, change Char to match labels
  if (!is.null(labels)) {
    char$Char <- labels[as.numeric(char$Char) + 1]
  }
  
  # return
  return(char)
}

###
# auxiliary function to get multipliers for discretized gamma
discretized_gamma_mult <- function(n_cat, alpha) {
  # bins for the categories
  a <- qgamma(seq(0, 1, 1 / n_cat)[-c(1, n_cat + 1)], alpha, alpha)
  
  # create multiplier vector
  r <- numeric(n_cat)
  
  # iterate through categories
  for (i in 1:n_cat) {
    # lower and upper limits for integration
    lower <- ifelse(i == 1, 0, a[i - 1])
    upper <- ifelse(i == n_cat, Inf, a[i])
    
    # get the expected value over that interval
    integral <- integrate(function(x) x * dgamma(x, alpha, alpha), 
                          lower, upper)$value
    
    # multiply by the number of categories, since each bin has prob 1/k
    r[i] <- n_cat * integral
  }
  
  # return multipliers
  return(r)
}

