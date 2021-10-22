#' Separate a paleobuddy simulation into monophyletic clades
#'
#' Separates a \code{sim} object into \code{sim} objects each with a mother
#' species and its descendants. If argument "S" is not used, it 
#' returns by default the list of \code{sim} objects
#' descended from each species with an \code{NA} parent in the original input 
#' (meaning species alive at the beginning of the simulation, \code{paleobuddy}
#' temrinology). Argument "S" allows user to input a vector of species whose 
#' daughters + itself will be the clade that makes the returned "sim" 
#' object (i.e., the subclade formed by "s" + daughters).
#' Returns for each clade a vector with the original 
#' identity of member species as well.
#'
#' @inheritParams make.phylo
#'
#' @param S A vector of species in \code{sim}. If not supplied, \code{S} will be
#' the starting species in the simulation (i.e. those for which the parent is
#' \code{NA} - please note that if only one species has \code{NA} as parent,
#' function will return the same object).
#'
#' @author Bruno Petrucci and Matheus Januario.
#'
#' @return A \code{list} object with (named) \code{sim} objects corresponding to 
#' the clades descended from species in \code{S}. For each clade, an extra vector 
#' \code{LIN} is included so the user can identify the order of species in the
#' returned \code{sim} object with the order of species in the original simulation.
#'
#' @examples
#'
#' # we will start with examples where S are the starting species
#' 
#' ###
#' # first, let us try a simulation with 3 clades,
#' set.seed(4)
#' sim <- bd.sim(n0 = 3, lambda = 0.1, mu = 0.1, tMax = 10, nFinal = c(20, Inf))
#' 
#' # using the functions
#' clades <- find.lineages(sim)
#' 
#' # set up for plotting side by syde
#' par(mfrow = c(1, length(clades)))
#' 
#' # for each clade
#' for (i in 1:length(clades)) {
#'   # change NA to 0 on the clade's TE
#'   clades[[i]]$sim$TE[clades[[i]]$sim$EXTANT] <- 0
#'   
#'   # if there is only one lineage in the clade
#'   if (length(clades[[i]]$sim$TE) < 2) {
#'     # placeholder plot
#'     plot(NA, xlim = c(-1, 1), ylim = c(-1, 1))
#'     text("simulation with \n just one lineage", x = 0, y = 0.5, cex = 2)
#'   }
#'   
#'   # if it is a proper phylogeny
#'   else {
#'     if (requireNamespace("ape", quietly = TRUE)) {
#'       plot <- ape::plot.phylo(
#'         make.phylo(clades[[i]]$sim),
#'         main = "red: extinction events \n blue: speciation events");
#'       ape::axisPhylo()
#'     }
#'     
#'     # checking speciation times
#'     for (j in 2:length(clades[[i]]$sim$TS)) {
#'       # the subtraction is just to adjust the wt with the plot scale
#'       lines(x = c(
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TS[j],
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TS[j]),
#'         y = c(plot$y.lim[1], plot$y.lim[2]), lwd = 2, col = "blue")
#'     }
#'     
#'     # checking extinction times:
#'     for (j in 1:length(sim$TE)) {
#'       # the subtraction is just to adjust the wt with the plot scale
#'       lines(x = c(
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TE[j],
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TE[j]),
#'         y = c(plot$y.lim[1], plot$y.lim[2]), lwd = 2, col = "red")
#'     }
#'   }
#' }
#' 
#' ###
#' # it works with any number of clades, of course
#' set.seed(8)
#' sim <- bd.sim(n0 = 5, lambda = 0.1, mu = 0.08, tMax = 10, nFinal = c(20, Inf))
#' 
#' # using the functions
#' clades <- find.lineages(sim)
#' 
#' # set up for plotting side by syde
#' par(mfrow = c(1, length(clades)))
#' 
#' # for each clade
#' for (i in 1:length(clades)) {
#'   # change NA to 0 on the clade's TE
#'   clades[[i]]$sim$TE[clades[[i]]$sim$EXTANT] <- 0
#'   
#'   # if there is only one lineage in the clade
#'   if (length(clades[[i]]$sim$TE) < 2) {
#'     # placeholder plot
#'     plot(NA, xlim = c(-1, 1), ylim = c(-1, 1))
#'     text("simulation with \n just one lineage", x = 0, y = 0.5, cex = 2)
#'   }
#'   
#'   # if it is a proper phylogeny
#'   else {
#'     if (requireNamespace("ape", quietly = TRUE)) {
#'       plot <- ape::plot.phylo(
#'         make.phylo(clades[[i]]$sim),
#'         main = "red: extinction events \n blue: speciation events");
#'       ape::axisPhylo()
#'     }
#'     
#'     # checking speciation times:
#'     for (j in 2:length(clades[[i]]$sim$TS)) {
#'       # the subtraction is just to adjust the wt with the plot scale
#'       lines(x = c(
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TS[j],
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TS[j]),
#'         y = c(plot$y.lim[1], plot$y.lim[2]), lwd = 2, col = "blue")
#'     }
#'     
#'     # checking extinction times:
#'     for (j in 1:length(sim$TE)) {
#'       # the subtraction is just to adjust the wt with the plot scale
#'       lines(x = c(
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TE[j],
#'         sort(clades[[i]]$sim$TS, decreasing = TRUE)[2] - 
#'            clades[[i]]$sim$TE[j]),
#'         y = c(plot$y.lim[1], plot$y.lim[2]), lwd = 2, col = "red")
#'     }
#'   }
#' }
#' 
#' ###
#' # including one clade
#' # To show the behavior of the function, we will not use the "S" argument.
#' As there is only one lineage whose parent is "NA", find.linages() will 
#' return the same phylogeny
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.08, tMax = 10, nFinal = c(5, Inf))
#' 
#' par(mfrow = c(1, 2))
#' 
#' # plotting sim and find.lineages(sim) - should be equal
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   ape::plot.phylo(make.phylo(sim), main="original")
#'   ape::axisPhylo()
#'   ape::plot.phylo(make.phylo(find.lineages(sim)[[1]]$sim), 
#'                   main="after find.lineages()")
#'   ape::axisPhylo()
#' }
#'
###
#' # now let us check that when S does not contain a starting species, we still
#' # get correct subsets of the simulation
#' set.seed(8)
#' sim <- bd.sim(1, 0.1, 0.05, 20, nFinal = c(5, Inf))
#' 
#' # making sure we have a couple of clades to explore
#' while ((length(which(sim$PAR == 1)) < 3) | (length(which(sim$PAR == 2)) < 3) |
#'        (length(which(sim$PAR == 3)) < 3)) {
#'   sim <- bd.sim(1, 0.1, 0.05, 20, nFinal = c(5, Inf))
#' }
#' 
#' par(mfrow=c(2,2))
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   # first we plot the clade started by 1
#'   ape::plot.phylo(make.phylo(sim), main="original")
#'   ape::axisPhylo()
#'   
#'   # this should look the same
#'   ape::plot.phylo(make.phylo(find.lineages(sim)[[1]]$sim), 
#'                   main="after find.lineages()")
#'   ape::axisPhylo()
#'   
#'   # and these should be part of the previous phylogenies
#'   ape::plot.phylo(make.phylo(find.lineages(sim, c(2, 3))$clade_2$sim),
#'                   main = "Daughters of sp 2")
#'   ape::axisPhylo()
#'   
#'   ape::plot.phylo(make.phylo(find.lineages(sim, c(2, 3))$clade_3$sim),
#'                   main = "Daughters of sp 3")
#'   ape::axisPhylo()
#' }
#' 
#' @name find.lineages
#' @rdname find.lineages
#' @export

find.lineages <- function(sim, S = NULL) {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # if S is null, the user wants to find the lineages with the simulation's
  # starting species as parents
  if (is.null(S)) {
    # by convention, species without parents in the output of the bd functions
    # have parents set to NA
    S = which(is.na(sim$PAR))
  }

  # create a final list
  final <- list()

  # find lineages for each species
  for (s in S) {
    # name the clade, and use the helper function below to find the species
    # descended from s for each s in S
    final[[paste0("clade_", s)]] = find.lineage(sim, s)
  }
  return(final)
}

###
# helper function for find.lineages

# does the exact same, but for one species

find.lineage <- function(sim, s) {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # if s is not on the simulation, we have a problem
  if (s > length(sim$TE)) {
    stop("This species is not on the simulation")
  }

  # lineage starts with a species
  lin <- c(s)

  # daughters of the first species of the lineage
  dau <- which(sim$PAR == s)

  # while species in the lineage have daughters
  while (length(dau) > 0) {
    # append the daughters to the lineage
    lin <- c(lin, dau)

    # find the daughters of the previous daughters
    dau <- which(sim$PAR %in% dau)
  }

  # make vectors for the clade
  TE <- sim$TE[lin]
  TS <- sim$TS[lin]
  PAR <- sim$PAR[lin]
  EXTANT <- sim$EXTANT[lin]

  # PAR here still follows the identifications on the original sim, so we need
  # to rename the species
  if (length(PAR) > 1) {
    # if the first species is not already unparented, it will be now
    PAR[1] = NA

    # first species of the clade (the one that generated the second) is 1
    PAR[PAR == PAR[2]] = 1

    # every other species follows the order in lin, to preserve the order
    # of TE and TS
    for (p in unique(PAR[PAR != 1 & !is.na(PAR)])) {
      PAR[PAR == p] = which(lin == p)
    }
  }

  # append it to a sim
  sim1 <- list(TE = TE, TS = TS, PAR = PAR, EXTANT = EXTANT)
  class(sim1) <- "sim"

  # note the inclusion of lin - this way, a user can tell which species in sim1
  # corresponded to which species in sim
  return(list(sim = sim1, LIN = lin))
}
