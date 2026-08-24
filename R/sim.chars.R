#' Discrete character simulation
#' 
#' Simulates a discrete character following a continuous-time Markov Chain
#' for all species in a \code{sim} and (optionally) a \code{sample} object.
#' The function simulates the character evolution for every species in 
#' \code{sim}, with each species taking as a starting state the state of its
#' parent at the time of speciation. For details on the character evolution
#' simulation algorithm, see \code{?sim.char.species}. Then, each character is
#' sampled for extant species and for each fossil sample, with options regarding
#' the state at which time should be assigned to a given species. Options are
#' also given regarding the labels of states, and rejection sampling for the
#' resulting data. Finally, the structure of parameters can accommodate
#' simulations up to a fixed number of characters, or simulations until we reach
#' a given number of characters sampled with given numbers of states.
#' 
#' @param sim A simulation object, usually returned by a birth-death simulation 
#' function like \code{bd.sim}.
#' 
#' @param sample A data frame with (at least) three columns: \code{Species} (the
#' species name, usually \code{tN.M} for the \code{M}th sample of the species 
#' that originated in \code{sim$TS[N]}), \code{Extant} (boolean indicating
#' whether the species is extant or not), and \code{SampT} (the fossil sampling
#' time of that species). Usually returned by a fossil sampling simulation
#' function like \code{sample.clade}. Default is \code{NULL}, in which case the
#' traits will only be sampled for the extant species in \code{sim}.
#' 
#' @param nTraits The number of traits to be simulated. It may be a number or
#' \code{NULL}. If it is a number, that number of traits will be simulated, with
#' \code{nStates} and \code{X0} expected to be either a single number or a
#' vector of length \code{nTraits} indicating the number of states for each
#' trait. If it is \code{NULL}, \code{nStates} is assumed to be a vector where
#' the nth element represents the desired number of traits with n states. This
#' is done because simulating a character with n states does not guarantee the
#' sampled data (at each specified time for each species) will lead to an n
#' state character. Setting \code{nTraits} and \code{nStates} as described will
#' guarantee \code{sum(nStates)} characters, each with the desired number of
#' states. Default is \code{1}.
#' 
#' @param nStates The number of states for each character. Can be either a
#' single number, a vector of numbers of length \code{nTraits}, or a vector of
#' numbers representing the number of characters to simulate with a given number
#' of states. See the \code{nTraits} description for details. Default is 
#' \code{2}.
#' 
#' @param X0 The initial state for each character. Can be either a single number
#' or a vector of numbers of length \code{nTraits}. If \code{nTraits} is
#' \code{NULL}, \code{X0} should be a single number representing the starting
#' state for all characters. Default is \code{0}.
#' 
#' @param Q The transition matrix for a given number of states. It is assumed
#' to be a function that takes one argument (the number of states) and returns
#' a square matrix with that many rows and columns. For different states 
#' \code{i} and \code{j}, the rate at which a species at state \code{i}
#' transitions to \code{j} is assumed to be \code{Q(n)[i + 1, j + 1]}, where 
#' \code{n} is the number of states. Default is a 2x2 matrix with symmetrical
#' transition rates equal to \code{2}. Note that while the standard is to set 
#' the diagonal values such that rows sum to \code{0}, the diagonal values are 
#' never used in our algorithm, and we generally set them to \code{0}.
#' 
#' @param species_clock The rate modified for each species in the simulation.
#' It can be a number (in which case it will be a global clock), or a vector of
#' length \code{length(sim$TS)}. Default is \code{1}, in which case the clock
#' will be defined by \code{Q}. That is the easiest way to set a strict clock.
#' For a relaxed clock, remember to set \code{Q} considering relative rates,
#' and \code{species_clock} will then be a multiplier for each species. Note
#' that since paleobuddy considers a budding speciation model, it is currently
#' impossible to have different rates at different branches belonging to the
#' same species. This might be implemented in the future.
#' 
#' @param var_min The percentage minimum variation for characters along samples.
#' After a character is simulated, the number of samples with each state up to
#' the number of states for that character is checked. If \code{nTraits} is 
#' \code{NULL}, it is assumed that only traits with exactly the given number of
#' states sampled will be taken, so a trait is rejected unless every state is
#' represented in at least \code{var_min*100} percent of species. If 
#' \code{nTraits} is a number, it is assumed that traits simulated with the 
#' given number of states, but that did not end up with that many states sampled
#' will be accepted. As such, then any trait with at least two states being
#' represented at this percentage or more will be accepted. For example, if
#' \code{var_min == 0} (the default), a 3-state trait will be accepted if either
#' all states are at all represented (if \code{nTraits == NULL}), or if at least
#' two states are (if \code{nTraits} is a number). In other words, only
#' invariant characters are rejected, but \code{nTraits} defines whether to
#' match the sampled and simulated number of states exactly.
#' 
#' @param labels The labels for the states on the output data. If \code{NULL}
#' (default), the states will be \code{0:(nStates - 1)}. If it is a vector, the
#' states will be changed to the labels, with state \code{0} being 
#' \code{labels[1]} and so forth. Note that currently it is assumed that one 
#' would use the same labels for all characters, so if multiple characters are
#' simulated with differing numbers of states, \code{labels} should be at least
#' as long as the number of states of the character with the most states. For
#' example, if a character with 2 states and one with 3 states are simulated,
#' and \code{labels} is \code{c("A", "B", "C")}, the first character will have
#' states \code{"A"} and \code{"B"}, while the second will have all three 
#' states.
#' 
#' @param sim_max The maximum number of simulations to run. If this number is
#' reached before we have simulated the requested number of traits (likely due
#' to \code{min_var} being too strict for the parameters given), the function
#' returns \code{NULL}. The default value is \code{10000}.
#' 
#' @param sample_fossil Which fossil to sample for each species' character. 
#' Acceptable values are \code{"first"}, \code{"last"}, and \code{"random"}.
#' If set to \code{"first"} or \code{"last"}, the respective occurrence time of
#' each species will be used to extract the character state for that species.
#' If set to \code{"random"}, each character will be taken from one of the
#' occurrence times, to be chosen by sampling all occurrence times at uniform
#' probability. This is necessary because the state of a character for a given
#' species may change throughout its lifetime. The default value is 
#' \code{"first"} for more predictable results, but a \code{"random"} value may
#' better emulate how morphological paleontology is carried in practice.
#' 
#' @param extant_sampling Whether to sample characters at the present for extant
#' species. Since extant samples often take precedence for the scoring of
#' morphological characters due to their completeness, one might want to sample
#' characters at the present for extant species regardless of the value of
#' \code{sample_fossil}. This argument is not necessary if 
#' \code{sample_fossil == "last"}. The default is \code{FALSE}, meaning extant
#' species will be sampled with the same procedure as extinct ones (defined by
#' \code{sample_fossil}), unless they have no fossil samples, in which case
#' characters must be sampled at the present.
#' 
#' @return A data frame with rows representing sampled species in the simulation
#' (i.e. extinct species with no fossil samples are excluded), and columns
#' representing each character.
#' 
#' @author Bruno do Rosario Petrucci
#' 
#' @references 
#' 
#' Jukes, T. H., Cantor, C. R. 1969. Evolution of protein molecules. New York:
#' Academic Press. 21-132.
#' 
#' Hasegawa, M., Kishino, H., Yano, T. 1985. Dating of the human-ape splitting 
#' by a molecular clock of mitochondrial DNA. Molecular Evolution. 
#' 22(2):160-174.
#' 
#' Lewis, P.O. 2001. A likelihood approach to estimating phylogeny from 
#' discrete morphological character data. Systematic Biology. 50(6):913-925.
#' 
#' @examples
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run simple simulations
#' sim <- bd.sim(1, 0.2, 0.05, 10)
#' sample <- suppressMessages(sample.clade(sim, 0.5, 10))
#' 
#' ###
#' # a simple example using all default parameters
#' 
#' # one character
#' chars <- sim.chars(sim, sample)
#' # note that this represents the states at which each species was at their
#' # first fossil occurrence times (if they were sampled)
#' 
#' ###
#' # we can make things more complicated
#' 
#' # let us simulate 10 traits
#' nTraits <- 10
#' 
#' # the first 5 will have 2 states, the last 5 will have 4
#' nStates <- c(rep(2, 5), rep(4, 5))
#' 
#' # we need to create a function for the state matrix
#' Q <- function(i) {
#'   res <- matrix(0.05, i, i)
#'   diag(res) <- 0
#'   return(res)
#' }
#' 
#' # run the simulation
#' chars <- sim.chars(sim, sample, nTraits, nStates, Q = Q)
#' # note that character 7 only had two states, even though we simulated 4
#' # this is fine, of course, but there is a way to force the simulation
#' # of characters with a number of observed states
#' 
#' # nStates is now a vector of how many characters of each state number we want
#' nStates <- c(0, 5, 0, 5)
#' 
#' # we simulate with nTraits = NULL
#' chars <- sim.chars(sim, sample, nTraits = NULL, nStates, Q = Q)
#' # now all characters have the observed number of states equal to the
#' # simulated number of states
#' # note that this might be a minor violation of the Mk model
#' 
#' ###
#' # species don't need to all have the same clock, of course
#' 
#' # return nStates to what it was previously
#' nStates <- c(rep(2, 5), rep(4, 5))
#' 
#' # we can modify Q to be a relative rate matrix
#' Q <- function(i) {
#'   res <- matrix(1 / (i - 1), i, i)
#'   diag(res) <- 0
#'   return(res)
#' }
#' 
#' # and draw species clocks from a lognormal
#' species_clock <- rlnorm(length(sim$TS), 0, 1)
#' 
#' # we will normalize it so that it has a mean of 1
#' species_clock <- species_clock / mean(species_clock)
#' 
#' # run the simulation
#' chars <- sim.chars(sim, sample, nTraits, nStates,
#'                    Q = Q, species_clock = species_clock)
#'                    
#' ###
#' # the labels parameter is useful when simulating DNA
#' 
#' # set the labels
#' labels <- c("A", "C", "G", "T")
#' 
#' # make a Q matrix corresponding to the Jukes-Cantor model
#' Q <- function(i) {
#'   res <- matrix(1/3, 4, 4)
#'   diag(res) <- 0
#'   return(res)
#' }
#' 
#' # run the simulation
#' chars <- sim.chars(sim, sample, nTraits = 10, nStates = 4, Q = Q,
#'                    labels = labels)
#' 
#' # you might then want to delete the data for the extinct species
#' chars[which(!sim$EXTANT), ] <- rep("?", 10)
#'                    
#' ###
#' # we can simulate under a more complicated molecular evolution model
#' # let us consider the HKY85 model (Hasegawa, Kishino and Yano 1985)
#' 
#' # let's assume a starting frequency slightly higher for A and C
#' pi_A <- pi_C <- 0.3
#' pi_G <- pi_T <- 0.2
#' 
#' # set starting frequency
#' X0 <- c(0, 0, 0, 1, 1, 2, 2, 2, 3, 3)
#' # in a real simulation study you would draw this, but we are just
#' # keeping it simple
#' 
#' # assume transitions happen twice as often as transverions
#' kappa <- 2
#' 
#' # set the Q matrix
#' Q <- function(i) return(matrix(c(0, pi_A, kappa * pi_A, pi_A,
#'                                  pi_C, 0, pi_C, kappa * pi_C,
#'                                  kappa * pi_G, pi_G, 0, pi_G,
#'                                  pi_T, kappa * pi_T, pi_T, 0), 4, 4))
#' # remember in R matrices are set column-wise
#' 
#' # run the simulation
#' chars <- sim.chars(sim, sample, nTraits = 10, nStates = 4, X0, Q, 
#'                    labels = labels)
#' chars[which(!sim$EXTANT), ] <- rep("?", 10)
#'                    
#' 
#' @name sim.chars
#' @rdname sim.chars
#' @export

sim.chars <- function(sim, sample = NULL, nTraits = 1, nStates = 2, X0 = 0, 
                      Q = function(i) return(matrix(c(0, 0.1, 0.1, 0), 2, 2)), 
                      species_clock = 1,
                      var_min = 0, 
                      labels = NULL, sim_max = 10000,
                      sample_fossil = "first", 
                      extant_sampling = FALSE) {
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
  
  # if sample is null, make it an empty data frame
  if (is.null(sample)) {
    sample <- data.frame(matrix(nrow = 0,
                                ncol = 3))
    colnames(sample) <- c("Species", "Extant", "SampT")
  }
  
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
                                 species_clock, labels,
                                 sample_fossil, extant_sampling)
          # allow for invariant character
          
          # increase sim_count
          sim_count <- sim_count + 1
          if (!is.null(sim_max) && sim_count > sim_max) return(NULL)
          
          # get states
          states <- if (is.null(labels)) 
            0:(i - 1) else 
              labels[1:i]
          
          # condition on minimum variation
          cond <- all(sapply(states, function(x) sum(char_n[, 2] == x)) / 
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
      
      # initialize condition to false
      cond <- FALSE
      
      # start sim counter--we don't want to try to run this for too long
      sim_count <- 1
      
      # until we have the number we want
      while (!cond) {
        # get this character 
        char_n <- sim.char.one(sim, sample, nStates[i], X0[i], Q_i,
                               species_clock, labels,
                               sample_fossil, extant_sampling)
        
        # increase sim_count
        sim_count <- sim_count + 1
        if (!is.null(sim_max) && sim_count > sim_max) return(NULL)
        
        # get states
        states <- if (is.null(labels)) 
          0:(nStates[i] - 1) else 
            labels[1:nStates[i]]
        
        # condition on minimum variation
        cond <- sum(sapply(states, 
                           function(x) sum(char_n[, 2] == x)) / 
                      nrow(char_n) > var_min) > 1
      }
      
      # cbind it to chars
      chars[, 1 + i] <- char_n$Char
    }
  }
  
  # make rownames species and delete first column
  rownames(chars) <- chars[, 1]
  chars <- chars[, -1, drop = FALSE]
  
  # make colnames just trait numbers
  colnames(chars) <- 1:ncol(chars)
  
  # return chars
  return(chars)
}

sim.char.one <- function(sim, sample, nStates, X0, Q, species_clock, 
                         labels = NULL, sample_fossil = "first",
                         extant_sampling = FALSE) {
  # check if Q has the correct number of dimensions
  if (nrow(Q) != nStates || nrow(Q) != nStates) {
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
  
  # if labels is not null, check that it is longer than the number of states
  if (!is.null(labels)) {
    if (length(labels) < nStates) {
      stop("There should be at least as many labels as states")
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
    
    # check if we want to sample at present
    if (sim$EXTANT[i] && (extant_sampling || nrow(sample_sp) == 0)) {
      char <- rbind(char, c(paste0("t", i), tail(traits[[i]]$value, 1)))
    } else if (nrow(sample_sp) > 0) {
      # if we don't want sampling at present,
      # check sample_fossil for which sample to get
      samp_to_get <- switch(sample_fossil,
                            first = 1,
                            last = length(sample_sp$SampT),
                            random = sample(1:length(sample_sp$SampT),
                                            1))
      # get time of a given fossil sample
      char_samp_time <- sample_sp$SampT[samp_to_get]
      
      # and sample state at that time (simulating a morphospecies)
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
  
  # remove unsampled species, if any
  unsampled <- which(is.na(char$Char))
  if (length(unsampled) > 0) char <- char[-unsampled, ]
  
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

