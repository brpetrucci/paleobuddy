#' Converts a paleobuddy simulation into a phylogeny
#'
#' \code{MakePhylo} generates a \code{phylo} object using a simulation from the
#' \code{BDSim} function. The phylogeny follows a "Hennigian" (sensu Ezard et
#' al 2011) format. If the simulation has only one lineage, the function
#' returns \code{NA} as there is no phylogeny for a simulation with only one
#' lineage.
#'
#' @param sim a simulation from the \code{BDSim} function.
#'
#' @author Function written by Matheus Januario. Reference: Ezard, T. H.,
#' Pearson, P. N., Aze, T., & Purvis, A. (2012). The meaning of birth and death
#' (in macroevolutionary birth-death models). Biology letters, 8(1), 139-142.
#'
#' @return A \code{phylo} object.
#'
#' @examples
#'
#' #Generating a phylogeny using constant rates
#' library(ape)
#' sim<-BDSim(N0 = 1, pp = 0.2, qq = 0.05, tmax = 10)
#' while(length(sim$TE) < 2){ #in case first simulation has only one species
#'   sim<-BDSim(N0 = 1, pp = 0.2, qq = 0.05, tmax = 10)
#' }
#' phy<-MakePhylo(sim)
#' plot.phylo(phy)
#' # we can also plot only the molecular phylogeny
#' plot.phylo(drop.fossil(phy))
#'
#' # this works for sim generated with any of the scenarios in \code{BDSim}, of course
#' sim<-BDSim(N0=1, pp=function(t) 0.12+0.01*t,qq=10, tmax=10, eshape=1.3)
#' while(length(sim$TE) < 2){ #in case first simulation has only one species
#'   sim<-BDSim(N0=1, pp=function(t) 0.12+0.01*t,qq=10, tmax=10, eshape=1.3)
#' }
#' phy<-MakePhylo(sim)
#' plot.phylo(phy)
#'
#' @name MakePhylo
#' @rdname MakePhylo
#' @export

MakePhylo<-function(sim){

  if(length(sim$TE)<2){
    message("There is no phylogeny for a simulation with only one lineage")
    return(NA)
  }

  all.dir.daugthers<-function(lin, x){
    #all.dir.daugthers returns the name of each direct daugther species
    #x = a simulation from PaleoBuddy
    #lin = a numeric specyfing the name of a lineage
    return(which(x$PAR ==lin))
  }


  cur.node<-length(sim$TE)+1 #current node
  edge<-matrix(nrow=1, ncol=2, data=c(cur.node, NA)) #edge matrix
  passed<-vector() #lineages which the function already put in the phylogeny
  i<-2 #current lineage
  lins<-c(1,2) #lineages which the function still has to solve (at least)
  jump<-0 #internal variable to help control the node function
  Nnode<-length(sim$TE)-1 #number of nodes in the phylogeny
  births_node<-rep(NA, times=length(sim$TE)) #vector storing the node correspondent to each birth
  births_node[2]<-cur.node
  counter<-0 #needed for debuging

  while(length(lins)>0){ #while some tip dont have a place in the phylogeny:

    (dau<-all.dir.daugthers(lin = i, x = sim))
    (dau<-dau[!(dau %in% passed)])

    if(is.numeric(dau) & length(dau)>0){ #if lineage has daughters

      if(jump==1){ #if a whole clade has very recently being put in the phylogeny
        (cur.node<-max(edge)+1)

        if(is.na(edge[nrow(edge),2])){
          edge[nrow(edge),2]<-cur.node
        } else{#Note: I think this never happens, but maintained here anyway just because caution
          edge<-rbind(edge,
                     matrix(nrow = 1, ncol = 2, data=c(prev.node, cur.node)))
        }
        births_node[dau[1]]<-cur.node
        jump<-0

      } else{ #if the current lineage is a non/monophyletic branch
        cur.node<-cur.node+1
        if(is.na(edge[nrow(edge),2])){
          edge[nrow(edge),2]<-cur.node
        } else{
          edge<-rbind(edge,
                     matrix(nrow = 1, ncol = 2, data=c(cur.node-1, cur.node)))
        }
        births_node[dau[1]]<-cur.node
      }

      edge<-rbind(edge,
                 matrix(nrow = 1, ncol = 2, data=c(cur.node, NA)))
      lins<-c(lins,dau[1])
      i<-lins[length(lins)]
    } #END if lineage has daughters


    if(is.numeric(dau) & length(dau)==0){ #if lineage has no daughters

      if(is.na(edge[nrow(edge),2])){
        edge[nrow(edge),2]<-i
      } else{
        edge<-rbind(edge,
                   matrix(nrow = 1, ncol = 2, data=c(max(edge[!(duplicated(edge[,1])|duplicated(edge[,1], fromLast=TRUE)),1]), i)))
      }
      passed<-c(passed, i)
      lins<-lins[-length(lins)]
      i<-lins[length(lins)]
    }  #END if lineage has no daughters

    if(sum(edge[,1] %in% cur.node)>1){ #this means that the function reached the end of the lineage of the cur.node
      #the warning here only "affects" a a condition which is never satisfied (jump when there is previous opened edge).
      #I supressed it because it is anoying
      suppressWarnings({prev.node= max(edge[!(duplicated(edge[,1])|duplicated(edge[,1], fromLast=TRUE)),1])})
      jump<-1
    }

    #registering bugs (if any)
    counter<-counter+1
    if(counter > 10*dim(edge)[1]){return("The function is lost and seems that it will not find a phylogeny. Please report this error and provide this simulation for debugging")}

  } #END while loop

  #calculating edge length:
  edge.length<-vector()
  for(i in 1:nrow(edge)){
    aux1<-edge[i,1]
    aux2<-edge[i,2]

    if(aux2<=length(sim$TE)){#if the branch is a tip:
      edge.length[i]<- sim$TS[which(births_node==aux1)]-sim$TE[aux2]
    } else{
      edge.length[i]<- sim$TS[which(births_node==aux1)]-sim$TS[which(births_node==aux2)]
    }

  }

  #Tyding all together to create the phylo object
  phy<-list(
    tip.label=paste0("t", 1:length(sim$TE)),
    edge=edge,
    edge.length=edge.length,
    Nnode=length(sim$TE)-1)
  class(phy)<-"phylo"

  return(phy)
}
