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
#' the order in which lineages should be drawn, from the bottom (i.e. the
#' first integer) to the upper side (#th integer, with # = number of lineages
#' in \code{sim}) of the figure. Default value of this parameter is "TS", so by
#' default species will be sorted by order of origination in the simulation.
#' 
#' @param showLabel A \code{logical} on whether to draw species labels (i.e. 
#' species 1 being t1, species 2 t2 etc.). Default is \code{TRUE}.
#' 
#' @param lwdLin The relative thickness/size of all elements (i.e., lines and 
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
#' 

draw.sim <- function (sim, fossils = NULL, sortBy = "TS", 
                      lwdLin = 4, showLabel = TRUE, ...) {
  # make NAs 0
  sim$TE[is.na(sim$TE)] <- 0
  
  # set default parameters for plot and par   
  # will only be used if user does not input parameters (...)
  args <- list(...)
  
  # suppress y axis
  if (!("yaxt" %in% names(args))) {
    oldpar <- par(no.readonly = TRUE)
    par(yaxt = "n")
  }
  
  # y axis name
  if (!("ylab" %in% names(args))) {
    formals(plot.default)$ylab <- "Simulated lineages"
  }
  
  # x axis name
  if (!("xlab" %in% names(args))) {
    formals(plot.default)$xlab <- "Time (Mya)"
  }
  
  # limits in x axis (to enclose all species' durations)
  if (!("xlim" %in% names(args))) {
    formals(plot.default)$xlim <- 
      c(max(sim$TS, na.rm = TRUE), min(sim$TE, na.rm = TRUE) - 1)
  }
  
  # limits in y axis
  if (!("ylim" %in% names(args))) {
    formals(plot.default)$ylim <- c(0, (length(sim$TE)+1))
  }
  
  # no frame
  if (!("frame.plot" %in% names(args))) {
    formals(plot.default)$frame.plot <- FALSE
  }
  
  # check inputs
  # sim object
  if (!is.sim(sim)) {
    stop("sim must be a valid sim object")
  }
  
  # fossil data frame
  if (!is.null(fossils)) {
    if (!((c("SampT") %in% colnames(fossils)) | 
          all(c("MaxT", "MinT") %in% colnames(fossils)))) {
      stop("fossils must contain either a SampT or both a MinT and MaxT
           columns. See ?draw.sim and ?sample.clade.")
    }
  }
  
  # check sortBy
  if (!class(sortBy) %in% c("character", "integer")) {
    stop("sortBy should be a character or a vector of integers.")
  }
  else if ((class(sortBy) == "integer") & !(all(1:length(sim$TE) %in% 
                                                unique(sortBy)))) {
    stop("sortBy must skip no lineage, and all lineages
         should have unique indices.")
  }
  
  # function for drawing transparency 
  makeTransparent <- function(someColor, alpha = 25) {
    # make color
    newColor <- col2rgb(someColor)
    
    # apply transparency to color
    apply(newColor, 2, function(curcoldata) {
      rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], 
          alpha = alpha, maxColorValue = 255)
    })
  }
  
  # function for drawing jittered rectangles
  jitter_foo <- function(x) {
    # create return
    jit <- vector()
    
    # jitter each rectangle
    while (length(jit) < length(x)) {
      # get a jitter amount
      r <- rnorm(1, sd = 2)
      
      # if too high, repeat
      while (r > 0.25 | r < (-0.25)) {
        r <- rnorm(1, sd = 2)
      }
      
      # append to result
      jit <- c(jit, r)
    }
    
    # return jittered rectangles
    return(x + jit)
  }
  
  # organize elements based on sortBy
  if (is.character(sortBy)) {
    # if it is PAR, need to be decreasing
    test <- c("PAR") %in% sortBy
    ord <- order(unlist(sim[sortBy]), decreasing = test)
  }
  else {
    # if it isn't a string, it's a vector of numbers
    ord <- sortBy
  }
  
  # copy of  sim object in the correct order
  sim_mod <- sim
  
  sim_mod$TE <- sim_mod$TE[ord]
  
  sim_mod$TS <- sim_mod$TS[ord]
  
  sim_mod$PAR <- sim_mod$PAR[ord]
  
  sim_mod$EXTANT <- sim_mod$EXTANT[ord]
  
  class(sim_mod) <- "sim"
  
  # open plot following user inputs + defaults above
  plot(NA, ...)
  
  # plot durations
  segments(x0 = sim_mod$TS, x1 = sim_mod$TE, y1 = 1:length(sim_mod$TE), 
           y0 = 1:length(sim_mod$TE), lwd = lwdLin, col = "black")
  
  # show species labels if user wants
  if (showLabel) {
    text(y = 1:length(sim_mod$TE), 
         x = sim_mod$TE - ((max(sim$TS) - min(sim$TE)) * 0.035), 
         labels = paste0("t", 
                         sprintf(paste0("%0", 
                                        ceiling(log(length(sim$TE), 10)), "d"), 
                                 ord)))
  }
  
  # establish references for plotting
  
  # find original species in the new sim object
  luca <- which(is.na(sim_mod$PAR))
  
  # find parents of each species in order
  aux_y <- unlist(lapply(sim$PAR, function(x) which(ord == x)[1]))
  
  # dashed lines between parent and daughter species
  segments(x0 = sim_mod$TS[-luca], x1 = sim_mod$TS[-luca], 
           y1 = (1:length(sim_mod$TE))[-luca], y0 = (aux_y[ord])[-luca], 
           lty = 2, lwd = lwdLin*0.25, col = "gray50")
  
  # add fossils
  if (!(is.null(fossils))) {
    # find number ID of each species' fossils
    ids <- as.numeric(gsub("t", "", fossils$Species))
    
    # if SampT is in the data frame, have true occurrence times
    if ("SampT" %in% colnames(fossils)) {
      # draw fossil occurrences as time points
      points(x = fossils$SampT, 
             y = unlist(lapply(ids, function(x) which(ord == x))), 
             col = "red", pch = 16, cex = lwdLin*0.25)
      
      # if MaxT and MinT are also in the columns, message
      if ("MaxT" %in% colnames(fossils) & "MinT" %in% 
          colnames(fossils)) {
        message("fossils contains both SampT and MaxT/MinT columns. Only
                true fossil occurrences will be drawn.")
      }
    }
    # if MaxT and MinT are in the data frame, have occurrence time ranges
    else if ("MaxT" %in% colnames(fossils) & "MinT" %in% 
             colnames(fossils)) {
      # jitter lines for drawing fossil ranges
      y_jittered <- jitter_foo(unlist(lapply(ids, 
                                             function(x) which(ord == x))))
      
      # draw fossil ranges a bit transparent
      segments(x1 = fossils$MaxT, x0 = fossils$MinT, y1 = y_jittered, 
               y0 = y_jittered, col = makeTransparent("red", 100), 
               lwd = lwdLin * 0.75)
    }
  }
 
#returning par to pre-function settings
if (!("yaxt" %in% names(args))){
  on.exit(par(oldpar))
 }
      
}
