#' Draw a sim object
#'
#' Draws species longevities for a paleobuddy simulation (a \code{sim} object -
#' see \code{?sim}) in the graphics window. Allows for the assignment of 
#' speciation and sampling events, and further customization.
#'
#' @inheritParams traits.summary
#'
#' @param sim A \code{sim} object, containing extinction times, speciation 
#' times, parent, and status information for each species in the simulation. 
#' See \code{?sim}.
#' 
#' @traits 
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
#' @param lineageColors Character vector giving the colors of all lineages, 
#' sorted by the original lineage order (the one in the \code{sim} object). 
#' Must have same length as the number of lineages in the \code{sim} object. 
#' If NULL (default value) all lineages are plotted as black. this parameter 
#' has no effect if \code{traits} is also provided.
#' 
#' @param tipLabels Character vector manually assigning the tip labels of all 
#' lineages, sorted by the original lineage order (the one in the \code{sim} 
#' object). Must have same length as the number of lineages in the \code{sim} 
#' object. If NULL (default value) all lineages are plotted as "t#", with "#" 
#' being the position of that lineage in the \code{sim} object.
#' 
#' @param traitID Numerical giving the trait which will be plotted. this 
#' parameter is only useful when multiple traits were simualted in the 
#' same \code{sim} object.
#' 
#' @param traitColors Character vector providing colors providing colors for 
#' the states of a given trait, so its length must equal or exceed the number 
#' of states. Default values provide 7 colors (and so they can plot up to 7 
#' states).
#' 
#' @param TraitLegendPlace Placement of state legend. Accepted values are 
#' \code{"topleft"} (default value), \code{"bottomleft"}, \code{"bottomright"}, 
#' and \code{"topright"}.
#' 
#' @param FossilsToDraw Character assigning if fossils will be represented by 
#' exact time placements (\code{"exact"}, efault value), by horizontal bars 
#' giving range information (\code{"ranges"}), or by both forms (\code{"all"}).
#' 
#' @param fossilRangeAlpha Numerical giving color transparency for fossil range 
#' representation. Integers between \code{0} and \code{255} are preferred, 
#' but any float between \code{0} and \code{1} is also accepted. Default value 
#' is \code{100}.
#' 
#' @param restoreOldPar Logical assigning if plot default values show be 
#' restored after function finalizes plotting. Deafult is \code{TRUE}, but 
#' users interesting in using plot additions (e.g. \code{abline()} to highlight 
#' a certain age) should assign this as \code{FALSE} to use the x and y values 
#' in the plot. If false, x-axis follows time, and y-axis follows the number of 
#' species plotted, with 1 being the bottom lineage, and the upper y-limit 
#' being the Nth lineage in the \code{sim}.
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
#' 
#' # we start drawing a simple simulation
#' #'
#' # maximum simulation time
#' tMax <- 10
#' #'
#' # set seed
#' set.seed(1)
#' #'
#' # run a simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#' 
#' # draw it
#' draw.sim(sim)
#' 
#' ###
#' # we can add fossils to the drawing
#' #'
#' # maximum simulation time
#' tMax <- 10
#' #'
#' # set seed
#' set.seed(1)
#' #'
#' # run a simulation
#' sim <- bd.sim(n0 = 1, lambda = 0.6, mu = 0.55, tMax = tMax, 
#'               nFinal = c(10,20)) 
#' #'
#' # set seed
#' set.seed(1)
#' #'
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
#' #'
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
#' draw.sim(sim, fossils = fossils, sortBy = "PAR",
#'          FossilsToDraw = "ranges", restoreOldPar = F)
#' 
#' # adding the bounds of the simulated bins
#' abline(v = bins, lty = 2, col = "blue", lwd = 0.5)
#' 
#' #alternatively, one can draw lineages varying colors and tip labels:
#' #(note how they are sorted)
#' draw.sim(sim, tipLabels = paste0("spp_", 1:10), 
#'          lineageColors = rep(c("red", "green", "blue"), times=5))
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
#' #'
#' # set seed
#' set.seed(1)  
#' #'
#' # simulate fossil sampling
#' fossils <- sample.clade(sim = sim, rho = 4, tMax = tMax, returnTrue = TRUE)
#' 
#' # draw it with random sorting (in pratice this could be a trait
#' # value, for instance)
#' draw.sim(sim, fossils = fossils, sortBy = sample(1:length(sim$TS)))
#' 
#' ###
#' # plotting trait dependent diversification:
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 20
#' 
#' # speciation, higher for state 1
#' lambda <- c(0.1, 0.2)
#' 
#' # extinction, lowest for state 0
#' mu <- c(0.01, 0.03)
#' 
#' # number of traits and states (2 binary traits)
#' nTraits <- 2
#' nStates <- 2
#' 
#' # initial value of both traits
#' X0 <- 0
#' 
#' # transition matrix, with symmetrical transition rates for trait 1,
#' # and asymmetrical (and higher) for trait 2
#' Q <- list(matrix(c(0, 0.1,
#'                    0.1, 0), ncol = 2, nrow = 2),
#'           matrix(c(0, 1,
#'                    0.5, 0), ncol = 2, nrow = 2))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nTraits = nTraits, 
#'                      nStates = nStates, X0 = X0, Q = Q, nFinal = c(2, 10))
#' 
#' # maybe we want to take a look at the traits of fossil records too
#' fossils <- sample.clade(sim$SIM, rho = 0.5, tMax = max(sim$SIM$TS), 
#'                         returnAll = T,bins = seq(0,20, by=1))
#' head(fossils)
#' 
#' draw.sim(sim$SIM, showLabel = T,  traits = sim$TRAITS, sortBy = "PAR",
#'          fossils = fossils, FossilsToDraw = "all")
#' 
#' 
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics points segments text
#' 
#' @name draw.sim
#' @rdname draw.sim
#' @export
#' 
draw.sim=function (sim, traits=NULL, fossils = NULL, lineageColors=NULL,
                   sortBy = "TS", lwdLin = 4, tipLabels=NULL, showLabel = TRUE,
                   traitID=1, traitColors=c("#a40000", "#16317d", "#007e2f", "#ffcd12", 
                                            "#b86092", "#721b3e", "#00b7a7"), TraitLegendPlace='topleft',
                   FossilsToDraw="exact", fossilRangeAlpha=100, restoreOldPar=TRUE, ...) {
  
  # set up to return par to pre-function settings
  oldPar <- par(no.readonly = TRUE) 
  
  if(restoreOldPar){
    on.exit(par(oldPar))  
  }
  
  # make NAs 0
  sim$TE[is.na(sim$TE)] <- 0
  
  # set default parameters for plot and par   
  # will only be used if user does not input parameters (...)
  args <- list(...)
  
  # suppress y axis
  if (!("yaxt" %in% names(args))) {
    #change par
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
    
    # warnings are not relevant here:
    suppressWarnings({
      min_time <- min(sim$TE, na.rm = TRUE)  
    })
    # if all species are extant, min_time must be adjusted:
    if(is.infinite(min_time)){
      min_time <- 0  
    }
    xmin <- min_time
    
    if(showLabel){
      xmin = 0 - (max(sim$TS, na.rm = TRUE)*0.035)
    }
    
    formals(plot.default)$xlim <- 
      c(max(sim$TS, na.rm = TRUE), xmin - 1)
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
    
    if(!(FossilsToDraw %in% c("all", "ranges", "exact"))){
      stop("\"FossilsToDraw \" must be equal to \"all\", \"ranges\", or to \"exact\"")
    }
    
    if (!((c("SampT") %in% colnames(fossils)) | 
          all(c("MaxT", "MinT") %in% colnames(fossils)))) {
      stop("fossils must contain either a SampT or both a MinT and MaxT
           columns. See ?draw.sim and ?sample.clade.")
    }
    
    if (!is.null(traits) & !(c("SampT") %in% colnames(fossils))){
      stop("Plotting fossils with traits require information on true times of fossilization. \n \n Try to create \"fossils\" using \"sample.clade(..., returnAll = T)\"") 
    }
  }
  
  
  # trait info:
  # getting traits:
  unique_trts <- vector()
  for(i in 1:length(traits)){
    unique_trts <- c(unique_trts, unique(traits[[i]][[traitID]]$value))
  }
  unique_trts <- unique(unique_trts)
  unique_trts
  
  if(length(traitColors)<length(unique_trts)){
    stop("Less colors than traits")
  }
  
  # checking length of inputted data:
  if(!is.null(traits)){
    unq_lengs <- length(unlist(unique(lapply(list(sim$TE, 
                                                  sim$TS, sim$PAR, sim$EXTANT, traits), length))))  
  }else{
    unq_lengs <- length(unlist(unique(lapply(list(sim$TE, 
                                                  sim$TS, sim$PAR, sim$EXTANT), length))))
  }
  
  if (unq_lengs != 1) {
    stop("Trait data does not match with birth-death simulation content. 
      Please check inputted sim and traits")
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
  
  #checking fossil input:
  if(!is.null(fossils)){
    if(!"SampT" %in% colnames(fossils) & FossilsToDraw != "ranges"){
      stop("No temporal information to draw requested fossil representation. 
           \n \n 
  Try to create \"fossils\" using \"sample.clade(..., returnAll = T)\".
  \n
  Alternatively, you might want to set FossilsToDraw = \"ranges\"")
    }
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
  
  # copy of sim object in the correct order
  sim_mod <- sim
  traits_mod <- traits
  
  sim_mod$TE <- sim_mod$TE[ord]
  
  sim_mod$TS <- sim_mod$TS[ord]
  
  sim_mod$PAR <- sim_mod$PAR[ord]
  
  sim_mod$EXTANT <- sim_mod$EXTANT[ord]
  
  traits_mod <-  traits_mod[ord]
  
  class(sim_mod) <- "sim"
  
  # open plot following user inputs + defaults above
  plot(NA, ...)
  
  # plot durations
  if(!is.null(traits)){
    for(i in 1:length(traits_mod)){
      segments(x0 = traits_mod[[i]][[traitID]]$max, 
               x1 = traits_mod[[i]][[traitID]]$min, 
               y1 = i, 
               y0 = i, 
               lwd = lwdLin, 
               # below we need the "+1" bc first trait is "zero"
               col = traitColors[traits_mod[[i]][[traitID]]$value+1])  
    }  
  }else{
    if(is.null(lineageColors)){
      lincols <- "black"
    }else{
      lincols <- lineageColors[ord]
    }
    
    segments(x0 = sim_mod$TS, x1 = sim_mod$TE, y1 = 1:length(sim_mod$TE), 
             y0 = 1:length(sim_mod$TE), lwd = lwdLin, col = lincols) 
  }
  
  
  
  
  # show species labels if user wants
  if (showLabel) {
    if(is.null(tipLabels)){
      tiplabels <- paste0("t", sprintf(paste0("%0", 
                                              ceiling(log(length(sim$TE), 10)), "d"), ord))
    }else{
      tiplabels <- tipLabels[ord]
    }
    
    labelPosition = sim_mod$TE
    labelPosition[is.na(labelPosition)]=0
    
    if(!is.null(traits)){
      textcolsaux = traits.summary(sim, traits = traits)[traitID]  
      textcols = traitColors[unlist(textcolsaux)+1]
      tip_traits <- unlist(lapply(traits, function(x) tail(x[[1]]$value, 1)))
      tip_cols = traitColors[tip_traits+1]
      tip_cols = tip_cols[ord]
    }else{
      tip_cols = rep("black", times=length(sim_mod$TE))
    }
    
    
    text(y = 1:length(sim_mod$TE), 
         x = labelPosition - ((max(sim$TS) - min_time ) * 0.035), 
         labels = tiplabels, col = tip_cols)
  }
  
  # establish references for plotting
  
  # find original species in the new sim object
  luca <- which(is.na(sim_mod$PAR))
  
  # find parents of each species in order
  aux_y <- unlist(lapply(sim$PAR, function(x) which(ord == x)[1]))
  
  # dashed lines between parent and daughter species
  if(!is.null(traits)){
    InheritedTraits = unlist(lapply(traits_mod, function(x) head(x[[traitID]]$value, 1)))
    vertLinCols = traitColors[InheritedTraits+1]
  }else{
    vertLinCols = "gray50"
  }
  
  
  segments(x0 = sim_mod$TS[-luca], x1 = sim_mod$TS[-luca], 
           y1 = (1:length(sim_mod$TE))[-luca], y0 = (aux_y[ord])[-luca], 
           lty = 2, lwd = lwdLin*0.25, col = vertLinCols)
  
  if(!is.null(traits)){
    legend(TraitLegendPlace, 
           legend= paste("State ", sort(unique(InheritedTraits)))
           , col=traitColors, lty=1, lwd=lwdLin)
    
  }
  
  # add fossils
  if (!(is.null(fossils))) {
    # find number ID of each species' fossils
    ids <- as.numeric(gsub("t", "", fossils$Species))
    
    if(!is.null(traits)){
      # stablishing traits:
      traitSummary <- traits.summary(sim, traits, 
                                     fossils = fossils, selection = "fossil")
      
      colFossils <- traitColors[unlist(traitSummary[traitID])+1]  
    }else{
      colFossils <- rep("red", times=length(sim$TE))
    }
    
    # if SampT is in the data frame, have true occurrence times
    if(FossilsToDraw=="ranges"){
      
      # jitter lines for drawing fossil ranges
      y_jittered <- jitter_foo(unlist(lapply(ids, 
                                             function(x) which(ord == x))))
      
      # draw fossil ranges a bit transparent
      segments(x1 = fossils$MaxT, x0 = fossils$MinT, y1 = y_jittered, 
               y0 = y_jittered, col = makeTransparent(someColor = colFossils,
                                                      alpha = fossilRangeAlpha), lwd = lwdLin * 0.75)
      
    }else{
      if(FossilsToDraw=="exact"){
        
        #raw fossil occurrences as time points
        points(x = fossils$SampT, 
               y = unlist(lapply(ids, function(x) which(ord == x))), 
               col = colFossils, pch = 16, cex = lwdLin*0.25)
        
      }else{
        
        # jitter lines for drawing fossil ranges
        y_jittered <- jitter_foo(unlist(lapply(ids, 
                                               function(x) which(ord == x))))
        
        # FIRST draw fossil ranges a bit transparent
        segments(x1 = fossils$MaxT, x0 = fossils$MinT, y1 = y_jittered, 
                 y0 = y_jittered, col = makeTransparent(someColor = colFossils,
                                                        alpha = fossilRangeAlpha), lwd = lwdLin * 0.75)
        
        # THEN draw fossil occurrences as time points
        points(x = fossils$SampT, 
               y = unlist(lapply(ids, function(x) which(ord == x))), 
               col = colFossils, pch = 16, cex = lwdLin*0.25)
        
      }
    }
  }
}
