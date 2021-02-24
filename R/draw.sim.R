#' Draw a sim object
#'
#' Draws a paleobuddy simulation (a \code{sim} object - please see 
#' \code{?sim}) in the graphics window. Allows for the assign of speciation or 
#' sampling events, and further customization.
#'
#' @param sim A \code{sim} object, containing extinction times, speciation 
#' times, parent, and status information for each species in the simulation. 
#' See \code{?sim}.
#' 
#' @param fossil_occ A \code{data.frame} containing the fossil occurrences of 
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
#' @return A draw of the simulation in the graphics window. If the 
#' \code{fossil_occ} data.frame is inputted, its format will dictate how 
#' fossil occurrences will be plotted. If \code{fossil_occ} have a column named
#'  "SampT" (i.e. the fossil sampling is known with exactitude), fossil 
#' occurrences are assigned as dots. If \code{fossil_occ} have two columns, 
#' each named as "MaxT" and "MinT" (i.e. the early and late stage bounds 
#' associated with each occurrence), fossil occurrences are 
#' represented as slightly jittered, semitransparent bars indicating the early 
#' and late bounds of each fossil occurrence.
#' 
#' @author Matheus Januario. 
#' 
#' 
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics points segments text
#' @name draw.sim
#' @rdname draw.sim
#' @export
#'
#' @examples
#' #### Example 1
#'
#'#maximum simulation time
#'tMax=10 
#'
#'#runing simulation (biological process)
#' sim=bd.sim(1, .6, .55, tMax = tMax, nFinal = c(10,20)) 
#'
#' #simulating data resulting from a fossilization process (in this case, 
#' # with a record of perfect resolution in time)
#' fdt=sample.clade(sim=sim, rho = 4, tMax = tMax, returnTrue = TRUE)
#' 
#' #ploting:
#' draw.sim(sim, fossil_occ = fdt, sort_by = "PAR")
#' 
#' #### Example 2
#' 
#' #maximum simulation time
#' tMax=10 
#' 
#' #runing simulation (biological process)
#' sim=bd.sim(1, .6, .55, tMax = tMax, nFinal = c(10,20)) 
#' 
#' #simulating data resulting from a fossilization process (in this case, 
#' # with a record that has limited time resolution - i.e. occurrences
#' # are binned in time)
#' 
#' #first lets create random bins:
#' bins=c(tMax, 0, runif(n = rpois(1, lambda = 6), min = 0, max = tMax))
#' #then lets simulate the fosisliation process:
#' fdt=sample.clade(sim=sim, rho = 2, tMax = tMax, returnTrue = FALSE, bins = bins)
#' 
#' #ploting (and this time sorting lineages by their parent lineage):
#' draw.sim(sim, fossil_occ = fdt, sort_by = "PAR")
#' #adding the bounds of the simulated bins:
#' abline(v=bins, lty=2, col="red", lwd=.5)
#' 
#' #### Example 3
#'  
#' #maximum simulation time
#' tMax=10 
#' 
#' #runing simulation (biological process)
#' sim=bd.sim(1, .6, .55, tMax = tMax, nFinal = c(10,20)) 
#'  
#' #ploting (and this time sorting lineages randomly):
#' fdt=sample.clade(sim=sim, rho = 4, tMax = tMax, returnTrue = TRUE)
#' 
#' #assignin a random order of lineages to be ploted:
#' draw.sim(sim, fossil_occ = fdt, sort_by = sample(1:length(sim$TE)))
#' 
draw.sim=function(sim, fossil_occ=NULL, sort_by="TS"){
  
  #checking input:
  if(!is.sim(sim)){
    stop("sim input is not a sim object")
  }
  
  if(!is.null(fossil_occ)){
    if(!((c("SampT") %in% colnames(fossil_occ)) | all(c("MaxT", "MinT") %in% colnames(fossil_occ)))){
      stop("fossil_occ is wrongly assigned. please see help(draw.sim)")
    } 
  }
  
  if(!class(sort_by) %in% c("character", "integer", "numeric")){
    stop("sort_by should be a character or a vector of integers. please see help(draw.sim)")
  }else if(class(sort_by) == "integer" & length(sort_by) != length(sim$TE)){
    stop("must have the same length than all elements in the sim object")
  }
  
  #aux function:
  makeTransparent<-function(someColor, alpha=25)
  {
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(
      red=curcoldata[1],
      green=curcoldata[2],
      blue=curcoldata[3],
      alpha=alpha, maxColorValue=255)})
  }
  
  jitter_foo=function(x){
    jit=vector()
      while(length(jit)<length(x)){
        r=1.2
        while(r>0.25 | r<(-0.25)){
          r=rnorm(1, sd = 2)
        }
        jit=c(jit, r)
      }
    return(x+jit)
  }
  
  #changking conventions
  sim$TE[is.na(sim$TE)]=0
  sim$TE[sim$TE<0]=0
  
  #reordering sim object
  if(is.character(sort_by)){
    if(sort_by %in% c("PAR")){
      test=TRUE
    }else{
      test=FALSE
    }
    
    ord=order(unlist(sim[sort_by]), decreasing = test)
  }else{
    ord=sort_by
  }
  
  sim_mod=sim
  sim_mod$TE=sim_mod$TE[ord]
  sim_mod$TS=sim_mod$TS[ord]
  sim_mod$PAR=sim_mod$PAR[ord]
  sim_mod$EXTANT=sim_mod$EXTANT[ord]
  
  #ploting
  plot(NA, xlim=c(max(sim$TS), min(sim$TE)-.6), ylim=c(1,max(ord)), yaxt="n", ylab="Simulated lineages", xlab="Time (Mya", frame.plot = F)
  
  #adding durations
  segments(x0 = sim_mod$TS, x1=sim_mod$TE, y1=1:length(sim_mod$TE),
           y0=1:length(sim_mod$TE), lwd=4, col="black")
  
  #adding labels
  text(y=1:length(sim_mod$TE), x=sim_mod$TE-((max(sim$TS) - min(sim$TE) )*.035), labels = paste0("t", sprintf(paste0("%0",round(length(sim$TE)/10, digits = 0),"d"), ord)))
  
  #adding budding events
  luca=which(is.na(sim_mod$PAR))
  aux_y=unlist(lapply(sim$PAR, function(x) which(ord==x)[1]))
  
  
  
  segments(x0 = sim_mod$TS[-luca], 
           x1 = sim_mod$TS[-luca], 
           y1 = (1:length(sim_mod$TE))[-luca], 
           y0 = (aux_y[ord])[-luca],
           lty = 2, lwd=1, col="gray50")
  
  
  #adding fossil occs:
  if(!(is.null(fossil_occ))){
    if("SampT" %in% colnames(fossil_occ)){
      ids=as.numeric(gsub('spp_', '', fossil_occ$Species))
      points(x=fossil_occ$SampT, 
             y=unlist(lapply(ids, function(x) which(ord==x))), col="red", pch=16)
    }else if("MaxT" %in% colnames(fossil_occ) & "MinT" %in% colnames(fossil_occ)){
      
      ids=as.numeric(gsub('spp_', '', fossil_occ$Species))
      
      y_jittered=jitter_foo(unlist(lapply(ids, function(x) which(ord==x))))
      
      segments(x1=fossil_occ$MaxT, x0 = fossil_occ$MinT, 
               y1=y_jittered,y0=y_jittered, 
               col=makeTransparent("red", 100), lwd=3)
      
    }else{
      stop("fossil_occ was wrongly inputed, please type help(draw.sim_mod)")
    } 
  }
  
}

