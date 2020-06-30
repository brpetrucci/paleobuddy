#' Wrapper for \code{BDSimConstant} and \code{BDSimGeneral}. Takes the most 
#' general set of rate-building vectors and functions, creates speciation and 
#' extinction rates using \code{MakeRate} and calls the appropriate function.
#'
#' \code{BDSim} takes an initial number of species, speciation and extinction
#' rate functions and a maximum time of simulation, together with multiple 
#' options to alter the rates, and calls \code{BDSimConstant} or 
#' \code{BDSimGeneral} with the rates.
#'
#' @param \code{N0} initial number of species, usually 1. Good param 
#' to tweak if one is observing a low sample size when testing.
#'
#' @param \code{pp} function to hold the speciation rate over time. It 
#' could be a function of time (to be an exponential rate or weibull scale), a
#' function of time and an environmental variable, or a vector of rates to be
#' accompanied by a vector of rate shifts \code{pshifts}. 
#'
#' @param \code{qq} similar as above, but for extinction rate.
#' 
#' Note: \code{pp} and \code{qq} must always be greater than 0
#'
#' @param \code{tmax} ending time of simulation. Any species still living 
#' after tmax is considered extant, and any species that would be generated 
#' after \code{tmax} is not born.
#'
#' @param \code{sshape} shape param for the Weibull distribution for 
#' age-dependent speciation. Default is 0, where \code{pp} will be considered a 
#' time-dependent exponential rate. For \code{sshape != NULL}, \code{pp} will
#' be considered a scale, and \code{rexp_var} will draw a Weibull distribution
#' instead.
#'
#' @param \code{eshape} similar as above, but for extinction rate.
#'
#' @param \code{env_pp} a matrix containing time points and values of an
#' enviromental variable, like temperature, for each time point. This will be
#' used to create a speciation rate, so \code{pp} must be a function of time 
#' and said variable.
#'
#' @param \code{env_qq} similar as above, but for extinction rate.
#'
#' @param \code{pshifts} vector of rate shifts. First element must be 
#' starting time for simulation (0). It must have the same length as \code{pp}.
#' E.g. \code{pp = c(0.1, 0.2, 0.1)}, \code{pshifts = c(0, 10, 20)} means p 
#' will be 0.1 from 0 to 10, 0.2 from 10 to 20, and 0.1 from 20 to \code{tmax}.
#' Note that using this  method for step-function rates is currently slower than using 
#' \code{ifelse}.
#' 
#' @param \code{qshifts} similar as above, but for extinction rate.
#'
#' @return the return list of either \code{BDSimConstant} or 
#' \code{BDSimGeneral}, which have the same elements, as follows
#' 
#' \describe{
#' \item{\code{TE}}{list of extinction times, with tmax + 0.01 as the time of 
#' extinction for extant species.}
#' 
#' \item{\code{TS}}{list of speciation times, with -0.01 as the time of 
#' speciation for species that started the simulation.}
#' 
#' \item{\code{PAR}}{list of parents. Species that started the simulation have 
#' NA, while species that were generated during the simulation have their 
#' parent's number. Species are numbered as they are born.}
#' 
#' \item{\code{EXTANT}}{list of booleans representing whether a species is 
#' extant.}}
#' 
#' @author written by Bruno do Rosario Petrucci.
#'
#' @examples
#'
#' since both \code{BDSimConstant} and \code{BDSimGeneral} have been tested in
#' their respective man pages, we will spend some time here giving examples of
#' possible combinations of diversification scenarios``

BDSim<-function(N0,pp,qq,tmax,sshape=NULL,eshape=NULL,env_pp=NULL,env_qq=NULL,
                pshifts=NULL,qshifts=NULL){
  # if we have ONLY numbers for pp and qq, it is constant
  if ((is.numeric(pp)&length(pp)==1)&
      (is.numeric(qq)&length(qq)==1&
       (is.null(c(sshape,eshape,env_pp,env_qq,pshifts,qshifts))))) {
    p<-pp
    q<-qq
    #' call BDSimConstant
    return(BDSimConstant(N0,p,q,tmax))
  }
  
  # else it is not constant
  else{
    # use MakeRate to create the rates we want
    p<-MakeRate(pp,env_pp,pshifts)
    q<-MakeRate(qq,env_qq,qshifts)
    
    # call BDSimGeneral
    return(BDSimGeneral(N0,p,q,tmax,sshape,eshape))
  }
}
