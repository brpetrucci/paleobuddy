#' Draw a sim object
#'
#' Draws species longevities for a paleobuddy simulation (a \code{sim} object -
#' see \code{?sim}) in the graphics window. Allows for the assignment of 
#' speciation and sampling events, and further customization.
#'
#' @param sim A \code{sim} object, containing extinction times, speciation 
#' times, parent, and status information for each species in the simulation. 
#' See \code{?sim}.
#' 
#' @param fossils A \code{data.frame} containing the fossil occurrences of 
#' each lineage, e.g. as returned by the \code{sample.clade} function. The
#' format of this argument will define the way fossils are drawn (see below).
#' 
#' @param sortBy A single character or integer vector indicating how lineages 
#' should be sorted in the plot. If it is a string (see example 3), it 
#' indicates which element in the \code{sim} object that should be used to sort
#' lineages in the plot. If it is a vector of integers, it directly specifies
#' the order in which lineages should be drawn (from the bottom 
#' (i.e. integer = 1)) to the upper side (integer = length of the \code{sim} 
#' elements) of the figure). Default value of this parameter is "TS", so by
#' default species will be sorted by order of origination in the simulation.
#' 
#' @param lwd_lin The relative thickness/size of all elements (i.e., lines and 
#' points in the plot. Default value is 4 (i.e. equal to \code{lwd = 4} for 
#' the black horizontal lines).
#' 
#' @param ... Further arguments to be passed to \code{plot}
#' 
#' @return A plot of the simulation in the graphics window. If the 
#' \code{fossils} data.frame is supplied, its format will dictate how fossil
#' occurrences will be plotted. If \code{fossils} has a \code{SampT} column
#' (i.e. the occurrence times are exact), fossil occurrences are assigned as 
#' dots. If \code{fossils} has columns \code{MaxT} and \code{MinT} (i.e. the 
#' early and late stage bounds associated with each occurrence), fossil 
#' occurrences are represented as slightly jittered, semitransparent bars 
#' indicating the early and late bounds of each fossil occurrence.
#' 
#' @author Matheus Januario
#'
#' @examples
#' 
#' ###
#' # we start drawing a simple simulation
#'
#' # maximum simulation time
#' tMax <- 10
#'
#' # set seed
#' set.seed(1)
#'
#' # run a simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#'               
#' # draw it
#' draw.sim(sim)
#' 
#' ###
#' # we can add fossils to the drawing
#'
#' # maximum simulation time
#' tMax <- 10
#'
#' # set seed
#' set.seed(1)
#'
#' # run a simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#'
#' # set seed
#' set.seed(1)
#'
#' # simulate data resulting from a fossilization process
#' # with exact occurrence times
#' fossils <- sample.clade(sim = sim, rho = 4, tMax = tMax, returnTrue = TRUE)
#' 
#' # draw it
#' draw.sim(sim, fossils = fossils)
#' 
#' # we can order the vertical drawing of species based on
#' # any element of sim
#' draw.sim(sim, fossils = fossils, sortBy = "PAR")
#' # here we cluster lineages with their daughters by
#' # sorting them by the "PAR" list of the sim object
#' 
#' draw.sim(sim, fossils = fossils, sortBy = "TE")
#' # here we sort lineages by their extinction times
#' 
#' ###
#' # try with fossil ranges
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run birth-death simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#' 
#' # simulate data resulting from a fossilization process
#' # with fossil occurrence time ranges
#' 
#' # set seed
#' set.seed(20)
#'
#' # create time bins randomly
#' bins <- c(tMax, 0, runif(n = rpois(1, lambda = 6), min = 0, max = tMax))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # simulate fossil sampling
#' fossils <- sample.clade(sim = sim, rho = 2, tMax = tMax, 
#'                         returnTrue = FALSE, bins = bins)
#' 
#' # draw it, sorting lineages by their parent
#' draw.sim(sim, fossils = fossils, sortBy = "PAR")
#' 
#' # adding the bounds of the simulated bins
#' abline(v = bins, lty = 2, col = "blue", lwd = 0.5)
#' 
#' ###
#' # we can control how to sort displayed species exactly
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run birth-death simulations
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#'
#' # set seed
#' set.seed(1)  
#'
#' # simulate fossil sampling
#' fossils <- sample.clade(sim = sim, rho = 4, tMax = tMax, returnTrue = TRUE)
#' 
#' # draw it with random sorting (in pratice this could be a trait
#' # value, for instance)
#' draw.sim(sim, fossils = fossils, sortBy = sample(1:length(sim$TS)))
#' 
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics points segments text
#' 
#' @name draw.sim
#' @rdname draw.sim
#' @export

draw.sim=function (sim, fossils = NULL, sortBy = "TS", lwd_lin=4, ...) {
  # make NAs 0
  sim$TE[is.na(sim$TE)] <- 0
  
  # set default parameters for plot and par
  args <- list(...)
  
  if (("yaxt" %in% names(args))) {
    par(yaxt = "n")
  }
  
  if (!("ylab" %in% names(args))) {
    formals(plot.default)$ylab <- "Simulated lineages"
  }
  
  if (!("xlab" %in% names(args))) {
    formals(plot.default)$xlab <- "Time (Mya)"
  }
  
  if (!("xlim" %in% names(args))) {
    formals(plot.default)$xlim <- 
      c(max(sim$TS, na.rm = TRUE), min(sim$TE, na.rm = TRUE) - 1)
  }
  
  if (!("ylim" %in% names(args))) {
    test <- sortBy %in% c("PAR")
    ord <- order(unlist(sim[sortBy]), decreasing = test)
    formals(plot.default)$ylim <- c(1, max(ord))
  }
  
  if (!("frame.plot" %in% names(args))) {
    formals(plot.default)$frame.plot <- FALSE
  }
  
  # check inputs

  if (!is.sim(sim)) {
    stop("sim must be a valid sim object")
  }
  
  if (!is.null(fossils)) {
    if (!((c("SampT") %in% colnames(fossils)) | all(c("MaxT", 
                                                      "MinT") %in% colnames(fossils)))) {
      stop("fossils must contain either a SampT or both a MinT and MaxT
           columns. See ?draw.sim and ?sample.clade.")
    }
  }
  
  if (!class(sortBy) %in% c("character", "integer", "numeric")) {
    stop("sortBy should be a character or a vector of integers.")
  }
  else if (class(sortBy) == "integer" & length(sortBy) != length(sim$TE)) {
    stop("sortBy must have the same length as elements of sim.")
  }
  else if ((class(sortBy) == "numeric") & !(all(1:length(sim$TE) %in% 
                                                unique(sortBy)))) {
    stop("sortBy must skip no lineage, and all lineages
         should have unique indices.")
  }
  
  makeTransparent <- function(someColor, alpha = 25) {
    newColor <- col2rgb(someColor)
    
    apply(newColor, 2, function(curcoldata) {
      rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], 
          alpha = alpha, maxColorValue = 255)
    })
  }
  
  jitter_foo <- function(x) {
    jit <- vector()
    
    while (length(jit) < length(x)) {
      r <- 1.2
      
      while (r > 0.25 | r < (-0.25)) {
        r <- rnorm(1, sd = 2)
      }
  
      jit <- c(jit, r)
    }
    
    return(x + jit)
  }
  
  if (is.character(sortBy)) {
    test <- sortBy %in% c("PAR")
    
    ord <- order(unlist(sim[sortBy]), decreasing = test)
  }
  else {
    ord <- sortBy
  }
  sim_mod <- sim
  
  sim_mod$TE <- sim_mod$TE[ord]
  
  sim_mod$TS <- sim_mod$TS[ord]
  
  sim_mod$PAR <- sim_mod$PAR[ord]
  
  sim_mod$EXTANT <- sim_mod$EXTANT[ord]
  
  plot(NA, ...)
  segments(x0 = sim_mod$TS, x1 = sim_mod$TE, y1 = 1:length(sim_mod$TE), 
           y0 = 1:length(sim_mod$TE), lwd = lwd_lin, col = "black")
  text(y = 1:length(sim_mod$TE), 
       x = sim_mod$TE - ((max(sim$TS) - min(sim$TE)) * 0.035), 
       labels = paste0("t", 
                       sprintf(paste0("%0", 
                                      round(log(length(sim$TE), 10), 
                                            digits = 0) + 1, "d"), ord)))
  
  luca <- which(is.na(sim_mod$PAR))
  
  aux_y <- unlist(lapply(sim$PAR, function(x) which(ord == x)[1]))
  
  segments(x0 = sim_mod$TS[-luca], x1 = sim_mod$TS[-luca], 
           y1 = (1:length(sim_mod$TE))[-luca], y0 = (aux_y[ord])[-luca], 
           lty = 2, lwd = lwd_lin*0.25, col = "gray50")
  
  if (!(is.null(fossils))) {
    ids <- as.numeric(gsub("t", "", fossils$Species))
    
    if ("SampT" %in% colnames(fossils)) {
      points(x = fossils$SampT, 
             y = unlist(lapply(ids, function(x) which(ord == x))), 
             col = "red", pch = 16, cex=lwd_lin*0.25)
    }
    
    else if ("MaxT" %in% colnames(fossils) & "MinT" %in% 
             colnames(fossils)) {
      y_jittered <- jitter_foo(unlist(lapply(ids, 
                                             function(x) which(ord == x))))
      
      segments(x1 = fossils$MaxT, x0 = fossils$MinT, y1 = y_jittered, 
               y0 = y_jittered, col = makeTransparent("red", 100), 
               lwd = lwd_lin * 0.75)
    }
  }
}
