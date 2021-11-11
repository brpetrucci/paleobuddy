#' Separate a paleobuddy simulation into monophyletic clades
#'
#' Separates a \code{sim} object into \code{sim} objects each with a mother
#' species and its descendants. If argument \code{S} is not used, it returns by
#' default the list of \code{sim} objects descended from each species with an
#' \code{NA} parent in the original input (meaning species alive at the 
#' beginning of the simulation). If a vector of numbers is supplied for 
#' \code{S}, the list of \code{sim} objects return will instead be descended 
#' from each species in \code{S}. Returns for each clade a vector with the 
#' original identity of member species as well.
#'
#' @inheritParams make.phylo
#'
#' @param S A vector of species in \code{sim}. If not supplied, \code{S} will be
#' the starting species in the simulation, i.e. those for which the parent is
#' \code{NA}. If only one species has \code{NA} as parent, there is only one
#' clade in the \code{sim} object, and therefore the function will return the 
#' input.
#'
#' @author Bruno do Rosario Petrucci and Matheus Januario.
#'
#' @return A \code{list} object with (named) \code{sim} objects corresponding 
#' to the clades descended from species in \code{S}. For each clade, an extra
#' vector \code{LIN} is included so the user can identify the order of species
#' in the returned \code{sim} objects with the order of species in the original
#' simulation.
#'
#' @examples
#' ###
#' # first, we run a simple simulation with one starting species
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run simulation with a minimum of 20 species
#' sim <- bd.sim(n0 = 3, lambda = 0.1, mu = 0.1, tMax = 10, 
#'               nFinal = c(20, Inf))
#'
#' # get a simulation object with the clade originating from species 2
#' clades <- find.lineages(sim, S = 2)
#' 
#' # now we can check to make sure the subclade was correctly separated
#' 
#' # change NA to 0 on the clade's TE
#' clades[[1]]$sim$TE[clades[[1]]$sim$EXTANT] <- 0
#' 
#' # plot the phylogeny
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   plot <- ape::plot.phylo(
#'     make.phylo(clades[[1]]$sim),
#'     main = "red: extinction events \n blue: speciation events");
#'   ape::axisPhylo()
#' }
#' 
#' # check speciation times
#' for (j in 2:length(clades[[1]]$sim$TS)) {
#'   # the subtraction is just to adjust the wt with the plot scale
#'   lines(x = c(
#'     sort(clades[[1]]$sim$TS, decreasing = TRUE)[2] -
#'        clades[[1]]$sim$TS[j],
#'     sort(clades[[1]]$sim$TS, decreasing = TRUE)[2] -
#'        clades[[1]]$sim$TS[j]),
#'     y = c(plot$y.lim[1], plot$y.lim[2]), lwd = 2, col = "blue")
#' }
#' 
#' # check extinction times:
#' for (j in 1:length(sim$TE)) {
#'   # the subtraction is just to adjust the wt with the plot scale
#'   lines(x = c(
#'     sort(clades[[1]]$sim$TS, decreasing = TRUE)[2] -
#'        clades[[1]]$sim$TE[j],
#'     sort(clades[[1]]$sim$TS, decreasing = TRUE)[2] -
#'        clades[[1]]$sim$TE[j]),
#'     y = c(plot$y.lim[1], plot$y.lim[2]), lwd = 2, col = "red")
#' }
#' 
#' ###
#' # now we try a simulation with 3 clades
#' 
#' # set seed
#' set.seed(4)
#' 
#' # run simulation
#' sim <- bd.sim(n0 = 3, lambda = 0.1, mu = 0.1, tMax = 10, 
#'               nFinal = c(20, Inf))
#' 
#' # get subclades descended from original species
#' clades <- find.lineages(sim)
#' 
#' # set up for plotting side by side
#' par(mfrow = c(1, length(clades)))
#' 
#' # for each clade
#' for (i in 1:length(clades)) {
#'   # change NA to 0 on the clade's TE
#'   clades[[i]]$sim$TE[clades[[i]]$sim$EXTANT] <- 0
#'   
#'   # if there is only one lineage in the clade, nothing happens
#'   if (length(clades[[i]]$sim$TE) < 2) {
#'     # placeholder plot
#'     plot(NA, xlim = c(-1, 1), ylim = c(-1, 1))
#'     text("simulation with \n just one lineage", x = 0, y = 0.5, cex = 2)
#'   }
#'   
#'   # else, plot phylogeny
#'   else {
#'     if (requireNamespace("ape", quietly = TRUE)) {
#'       plot <- ape::plot.phylo(
#'         make.phylo(clades[[i]]$sim),
#'         main = "red: extinction events \n blue: speciation events");
#'       ape::axisPhylo()
#'     }
#'     
#'     # check speciation times
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
#'     # check extinction times:
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
#' # we can also have an example with more non-starting species in S
#' 
#' # set seed
#' set.seed(3)
#' 
#' # run simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.1, tMax = 10, 
#'               nFinal = c(10, Inf))
#' 
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   # set up for plotting side by side
#'   par(mfrow = c(1, 2))
#'   
#'   # first we plot the clade started by 1
#'   ape::plot.phylo(make.phylo(sim), main = "original")
#'   ape::axisPhylo()
#'   
#'   # this should look the same
#'   ape::plot.phylo(make.phylo(find.lineages(sim)[[1]]$sim), 
#'                   main="after find.lineages()")
#'   ape::axisPhylo()
#'   
#'   # get sublcades descended from the second and third species
#'   clades <- find.lineages(sim, c(2,3))
#'   
#'   # and these should be part of the previous phylogenies
#'   ape::plot.phylo(make.phylo(clades$clade_2$sim),
#'                   main = "Daughters of sp 2")
#'   ape::axisPhylo()
#'   
#'   ape::plot.phylo(make.phylo(clades$clade_3$sim),
#'                   main = "Daughters of sp 3")
#'   ape::axisPhylo()
#' }
#' ###
#' # if there is only one clade and we use the default for 
#' # S, we get back the original simulation object
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.1, mu = 0.08, tMax = 10, 
#'               nFinal = c(5, Inf))
#'               
#' # set up for plotting side by side
#' par(mfrow = c(1, 2))
#' 
#' # plotting sim and find.lineages(sim) - should be equal
#' if (requireNamespace("ape", quietly = TRUE)) {
#' 
#'   ape::plot.phylo(make.phylo(sim), main="original")
#'   ape::axisPhylo()
#'   ape::plot.phylo(make.phylo(find.lineages(sim)[[1]]$sim), 
#'                   main="after find.lineages()")
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
