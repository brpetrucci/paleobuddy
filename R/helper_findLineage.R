#' Helper function for findLineages
#'
#' \code{helper_findLineage} takes a simulation object, usually from \code{BDSim},
#' and a species number, and returns the simulation object that has that
#' species as a common ancestor.
#'
#' @param sim a simulation from the \code{BDSim} function. The
#' function accept simulations with any number of starting species.
#'
#' @param s a species number. The function will return a sim object whose common
#' ancestor is the species s, and s will be renumbered to species 1 in the result.
#'
#' @author written by Bruno Petrucci.
#'
#' @return a \code{sim} object with \code{s} as the mother species, including
#' all descendants of \code{s}.
#'
#' Note: \code{helper_findLineage(sim, s)} will be identical to
#' \code{findLineages(sim, c(s))}, so we omit examples here.
#'
#' @export

helper_findLineage <- function(sim, s) {
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
  if(length(PAR)>1){
    # if the first species is not already unparented, it will be now
    PAR[1] = NA

    # first species of the clade (the one that generated the second) is 1
    PAR[PAR==PAR[2]] = 1

    # every other species follows the order in lin, to preserve the order
    # of TE and TS
    for (p in unique(PAR[PAR != 1 & !is.na(PAR)])) {
      PAR[PAR==p] = which(lin==p)
    }
  }

  # append it to a sim
  sim1 <- list(TE=TE, TS=TS, PAR=PAR,
               EXTANT=EXTANT)

  return(sim1)
}
