#' Draw a sim object
#'
#' Draws a paleobuddy simulation (a \code{sim} object - please see 
#' \code{?sim}) in the graphics window. Allows for the assignment of speciation
#' or sampling events, and further customization.
#'
#' @param sim A \code{sim} object, containing extinction times, speciation 
#' times, parent, and status information for each species in the simulation. 
#' See \code{?sim}.
#' 
#' @param fossils A \code{data.frame} containing the fossil occurrences of 
#' each lineage, as returned by the \code{sample.clade} function.
#' 
#' @param sort_by A single character or integer vector indicating how lineages 
#' should be sorted in the plot. If inputted as character (see example 3), it 
#' indicates which element in the \code{sim} object that should be used to sort 
#' lineages in the plot. If inputted as a vector of integers, it directly 
#' specifies the order in which lineages should be drawn (from the bottom 
#' (i.e. integer = 1)) to the upper side (integer = length of the \code{sim} 
#' elements) of the figure). Default value of this parameter is "TS".
#' 
#' @return A plot of the simulation in the graphics window. If the 
#' \code{fossils} data.frame is inputted, its format will dictate how fossil
#' occurrences will be plotted. If \code{fossils} has a \code{SampT} column
#' (i.e. the occurrence times are exact), fossil occurrences are assigned as dots. 
#' If \code{fossils} has columns \code{MaxT} and \code{MinT} (i.e. the early and 
#' late stage bounds associated with each occurrence), fossil occurrences are 
#' represented as slightly jittered, semitransparent bars indicating the early 
#' and late bounds of each fossil occurrence.
#' 
#' @author Matheus Januario
#'
#' @examples
#' 
#' ###
#' # start simple
#' 
#' # maximum simulation time
#' tMax <- 10 
#'
#' # run birth-death simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#'
#' # simulate data resulting from a fossilization process
#' # with exact occurrence times
#' fdt <- sample.clade(sim = sim, rho = 4, tMax = tMax, returnTrue = TRUE)
#' 
#' # draw it
#' draw.sim(sim, fossils = fdt, sort_by = "PAR")
#' 
#' ###
#' # try with fossil ranges
#' 
#' # maximum simulation time
#' tMax <- 10 
#' 
#' # run birth-death simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#' 
#' # simulate data resulting from a fossilization process
#' # with fossil occurrence time ranges
#' 
#' # create time bins randomly
#' bins <- c(tMax, 0, runif(n = rpois(1, lambda = 6), min = 0, max = tMax))
#' 
#' # simulate fossil sampling
#' fdt <- sample.clade(sim = sim, rho = 2, tMax = tMax, 
#'                     returnTrue = FALSE, bins = bins)
#' 
#' # draw it, sorting lineages by their parent
#' draw.sim(sim, fossils = fdt, sort_by = "PAR")
#' 
#' # adding the bounds of the simulated bins
#' abline(v = bins, lty = 2, col = "red", lwd = 0.5)
#' 
#' ###
#' # we can control how to sort displayed species exactly
#' 
#' # maximum simulation time
#' tMax <- 10 
#' 
#' # run birth-death simulations
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#'  
#' # simulate fossil sampling
#' fdt <- sample.clade(sim = sim, rho = 4, tMax = tMax, returnTrue = TRUE)
#'  
#' # draw it with random sorting
#' draw.sim(sim, fossils = fdt, sort_by = sample(1:length(sim$TE)))
#' 
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics points segments text
#' 
#' @name draw.sim
#' @rdname draw.sim
#' @export
#' 

draw.sim <- function(sim, fossils = NULL, sort_by = "TS"){
  # some error checking
  
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("sim must be a valid sim object")
  }
  
  # check that fossils has a SampT or both a MaxT and MinT columns
  if (!is.null(fossils)) {
    if (!((c("SampT") %in% colnames(fossils)) | 
          all(c("MaxT", "MinT") %in% colnames(fossils)))) {
      stop("fossils must contain either a SampT or both a MinT and MaxT
           columns. See ?draw.sim and ?sample.clade.")
    } 
  }
  
  # check that sort_by has an accepted type
  if (!class(sort_by) %in% c("character", "integer", "numeric")) {
    stop("sort_by should be a character or a vector of integers.")
  } else if (class(sort_by) == "integer" & length(sort_by) != length(sim$TE)) {
    
    # if it is a vector of integers with the wrong length, error
    stop("sort_by must have the same length as elements of sim.")
  } else if((class(sort_by) == "numeric") & 
            !(all(1:length(sim$TE) %in% unique(sort_by)))) {
    
    # if it s is numeric, it must include every lineage
    stop("sort_by must skip no lineage, and all lineages should have unique indices.")
  }
  
  # aux functions
  makeTransparent <- function(someColor, alpha = 25) {
    newColor <- col2rgb(someColor)
    apply(newColor, 2, function(curcoldata) {
      rgb(
      red = curcoldata[1],
      green = curcoldata[2],
      blue = curcoldata[3],
      alpha = alpha, maxColorValue = 255)
      })
  }
  
  jitter_foo <- function(x){
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
  
  # changing conventions
  sim$TE[is.na(sim$TE)] <- 0
  
  # reordering sim object
  if (is.character(sort_by)) {
    if (sort_by %in% c("PAR")) {
      test <- TRUE
    } else {
      test <- FALSE
    }
    
    ord <- order(unlist(sim[sort_by]), decreasing = test)
  } else {
    ord <- sort_by
  }
  
  sim_mod <- sim
  sim_mod$TE <- sim_mod$TE[ord]
  sim_mod$TS <- sim_mod$TS[ord]
  sim_mod$PAR <- sim_mod$PAR[ord]
  sim_mod$EXTANT <- sim_mod$EXTANT[ord]
  
  # plotting
  plot(NA, xlim = c(max(sim$TS), min(sim$TE) - 0.6), ylim = c(1, max(ord)), 
       yaxt = "n", ylab = "Simulated lineages", xlab = "Time (Mya)",
       frame.plot = FALSE)
  
  # adding durations
  segments(x0 = sim_mod$TS, x1 = sim_mod$TE, y1 = 1:length(sim_mod$TE),
           y0 = 1:length(sim_mod$TE), lwd = 4, col = "black")
  
  # adding labels
  text(y = 1:length(sim_mod$TE), 
       x = sim_mod$TE - ((max(sim$TS) - min(sim$TE)) * 0.035), 
       labels = paste0("t", 
                       sprintf(paste0("%0", 
                                      round(length(sim$TE) / 10, 
                                            digits = 0), "d"), ord)))
  
  # adding budding events
  luca <- which(is.na(sim_mod$PAR))
  aux_y <- unlist(lapply(sim$PAR, function(x) which(ord == x)[1]))
  
  segments(x0 = sim_mod$TS[-luca], 
           x1 = sim_mod$TS[-luca], 
           y1 = (1:length(sim_mod$TE))[-luca], 
           y0 = (aux_y[ord])[-luca],
           lty = 2, lwd=1, col="gray50")
  
  # adding fossils
  if (!(is.null(fossils))) {
    if ("SampT" %in% colnames(fossils)) {
      ids <- as.numeric(gsub('spp_', '', fossils$Species))
      points(x = fossils$SampT, 
             y =  unlist(lapply(ids, function(x) 
               which(ord == x))), col = "red", pch = 16)
    } else if("MaxT" %in% colnames(fossils) & 
              "MinT" %in% colnames(fossils)) {
      ids = as.numeric(gsub('spp_', '', fossils$Species))
      
      y_jittered <- jitter_foo(unlist(lapply(ids, function(x) 
                                                  which(ord == x))))
      
      segments(x1 = fossils$MaxT, x0 = fossils$MinT, 
               y1 = y_jittered, y0 = y_jittered, 
               col = makeTransparent("red", 100), lwd = 3)
      
    }
  }
  
}

