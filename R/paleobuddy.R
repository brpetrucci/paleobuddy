#' paleobuddy: Simulating diversification dynamics
#'
#' \code{paleobuddy} provides users with flexible scenarios for species 
#' birth-death simulation and fossil record generation, besides the possibility
#' of generating phylogenetic trees from the same underlying process. 
#' 
#' @section Birth-death simulation:
#' Users have access to a large array of scenarios to use and combine for species
#' birth-death simulation. The function \code{bd.sim} allows for constant rates,
#' rates varying as a function of time, or time and/or an environmental variable,
#' besides age-dependent rates with the use of a shape parameter. Extinction and 
#' speciation rates can be supplied independently, so that one can combine any
#' types of scenarios for birth and death rates. See \code{?bd.sim},
#' \code{?bd.sim.constant} and \code{?bd.sim.general} for more information.
#' 
#' @section Fossil record simulation:
#' The package provides users with a similarly diverse array of scenarios for 
#' preservation rates in generating fossil records from birth-death simulations.
#' The function \code{sample.clade} accepts constant and time-varying rates. 
#' Users might also supply a constant rate and an age-dependent preservation
#' model describing the distribution of fossil occurrences over a species
#' duration. We supply the function \code{find.lineages} to allow for birth-death
#' simulations that start with multiple species by separating those in
#' monophyletic clades so that one can generate fossil records and phylogenies
#' (see below) for clades with a specific mother species. See 
#' \code{?sample.clade}, \code{?sample}, \code{?sample.adpp} and
#' \code{?find.lineages} for more information.
#' 
#' @section Phylogeny generation:
#' We believe it is imperative to be able to generate fossil records and
#' phylogenetic trees from the same underlying process, so the package provides
#' \code{make.phylo}, a function that takes a simulation object of the form 
#' returned by \code{bd.sim} and generates a \code{phylo} object from the APE
#' package. One can then use functions such as \code{ape::plot.phylo} and
#' \code{ape::drop.fossil} to plot the phylogeny or analyze the molecular
#' phylogeny. Since APE is not required for any function in the package, it is
#' a suggested but not imported package. Note that, as above, the function
#' \code{find.lineages} allows users to separate clades with mother species of 
#' choice, the results of which can be passed to \code{make.phylo} to generate 
#' separate phylogenies for each clade. See \code{?make.phylo} and 
#' \code{?find.lineages} for more information.
#' 
#' @section Utility functions:
#' The package makes use of a few helper functions for simulating and testing
#' that we provide the user for completion. \code{rexp.var} aims to emulate the
#' behavior of \code{rexp}, with the possibility of supplying a time-varying
#' rate. The function also allows for a shape parameter, in which case the times
#' drawn will be distributed as a Weibull, for age-dependent rates. 
#' \code{var.rate.div} calculates the expected diversity of a birth-death process
#' with varying rates for any time period, and is useful when testing birth-death
#' simulation functions. Finally, \code{binner} simply returns the number of 
#' fossil occurrences in each time bin for an occurrence and a bins vector
#' supplied by the user. This is mostly for use in the \code{sample.clade}
#' function. See \code{?rexp.var}, \code{?var.rate.div} and \code{?binner} for
#' more information.
#' 
#' @author 
#' Bruno do Rosario Petrucci, Matheus Januario and Tiago B. Quental
#' 
#' Maintainer: Bruno do Rosario Petrucci <petrucci@iastate.edu>
#' 
#' @examples 
#'
#' # speciation rate
#' p <- function(t) {
#'   0.1 + 0.04*t
#' }
#' 
#' # extinction rate
#' q <- 0.08
#' 
#' # these are pretty simple scenarios, of course
#' # check the examples in ?bd.sim for a more comprehensive review
#' 
#' # diversification
#' d <- function(t) {
#'   p(t) - q
#' }
#' 
#' # calculate how many species we expect over 10 million years
#' # note we are starting with 3 species
#' div <- var.rate.div(ff = d, n0 = 3, t = seq(0, 10, 0.1))
#' 
#' # plot it
#' plot(seq(0, 10, 0.1), div, type = 'l', main = "Expected diversity",
#'      xlab = "Time (My)", ylab = "Species")
#' 
#' set.seed(3)
#' 
#' # around 28 species by the end, seems pretty good
#' # run the simulation
#' sim <- bd.sim(n0 = 3, pp = p, qq = q, tMax = 10, nFinal = c(20, Inf))
#' # nFinal controls the final number of species
#' # here we are telling bd.sim to throw away simulations with less than 20
#' 
#' set.seed(1)
#' 
#' # from sim, we can create fossil records for each species
#' samp <- sample.clade(sim = sim, rr = 0.75, tMax = 10,
#'                      bins = seq(10, 0, -1))
#' # note 15 out of the 29 species did not leave a fossil - we can in this way
#' # simulate the incompleteness of the fossil record
#' 
#' # take a look at the resulting data frame
#' samp
#' 
#' # we can separate it in three monophyletic clades
#' clades <- find.lineages(sim = sim)
#' # note if we wanted to check for clades originated from specific species we
#' # need only pass S as a list with those species to find.lineages
#' 
#' # get a phylogeny for whichever one has a lot of species
#' for (c in clades) {
#'   if (length(c$TE) > 10) {
#'     phy <- make.phylo(c)
#'     break
#'   }
#' }
#' 
#' # take a look at the phylogeny
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   par(mfrow = c(1, 2))
#'   ape::plot.phylo(phy)
#'   
#'   # we can also plot the molecular phylogeny
#'   ape::plot.phylo(ape::drop.fossil(phy))
#' }
#' 
#' @docType package
#' @name paleobuddy
NULL