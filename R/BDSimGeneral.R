#' Returns information of a simulated clade for general speciation and
#' extinction rates
#'
#' \code{BDSimGeneral} takes an initial number of species, speciation and
#' extinction rates (either functions of time or of time and some
#' environmental variable), a maximum simulation time and possibly a shape for
#' age-dependent speciation and/or extinction.
#'
#' @param N0 initial number of species, usually 1. Good param to
#' tweak if one is observing a low sample size when testing.
#'
#' @param pp function to hold the speciation rate over time.
#' \code{BDSim} supplies this function with a \code{pp} ready to be used, so
#' that the only other information \code{BDSimGeneral} needs is a shape in case
#' the rate is to be age-dependent.
#'
#' @param qq similar to above, but for extinction rate.
#'
#' @param tmax ending time of simulation. Any species still living
#' after \code{tmax} is considered extant, and any species that would be
#' generated after \code{tmax} is not born.
#'
#' @param pshape shape param for the Weibull distribution for
#' age-dependent speciation. Default is 0, where \code{pp} will be considered a
#' time-dependent exponential rate. For \code{pshape != NULL}, \code{pp} will
#' be considered a scale, and \code{rexp_var} will draw a Weibull distribution
#' instead.
#'
#' @param qshape similar as above, but for extinction rate.
#'
#' @param fast when \code{TRUE}, sets \code{rexp_var} to throw away waiting times
#' higher than the maximum simulation time. Should be \code{FALSE} for unbiased
#' testing. User might also se it to \code{FALSE} for more accurate waiting times.
#'
#' @param
#'
#' @return a list of vectors, as follows
#'
#' \describe{
#' \item{\code{TE}}{list of extinction times, with -0.01 as the time of
#' extinction for extant species.}
#'
#' \item{\code{TS}}{list of speciation times, with tmax+0.01 as the time of
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
#' # REMEMBER TO INVERT TE AND TS IN SIMMEAN()
#'
#' @name BDSimGeneral
#' @rdname BDSimGeneral
#' @export

BDSimGeneral<-function(N0,pp,qq,tmax,pshape=NULL,qshape=NULL,fast=TRUE){
  # create vectors to hold times of speciation, extinction, parents and status
  TS<-rep(-0.01,N0)
  TE<-rep(NA,N0)
  Parent<-rep(NA,N0)
  is.extant<-rep(TRUE,N0)

  # initialize species count
  Scount<-1

  # while we have species to be analyzed still
  while (length(TE)>=Scount){

    # get the time of speciation, or 0 if the species
    # was there at the beginning
    tNow<-ifelse(TS[Scount]<0,0,TS[Scount])

    # find the waiting time using rexp_var - note that in rexp_var we only
    # count t from tNow (to consider the rates as functions), so that
    # now we need to subtract tNow
    WaitTimeS<-ifelse(pp(tNow)>0,
                      rexp_var(1,pp,tNow,tmax,pshape,
                               ifelse(TS[Scount]<0,0,TS[Scount]), fast),Inf)
    WaitTimeE<-ifelse(qq(tNow)>0,
                      rexp_var(1,qq,tNow,tmax,qshape,
                               ifelse(TS[Scount]<0,0,TS[Scount]), fast),Inf)

    # if the time of extinction is after the end of the simulation, make it tmax
    tExp<-min(tNow+WaitTimeE, tmax)

    # while there are fast enough speciations before the species goes extinct,
    while ((tNow+WaitTimeS)<=tExp){

      # advance to the time of speciation
      tNow<-tNow+WaitTimeS

      # add new times to the vectors
      TS<-c(TS,tNow)
      TE<-c(TE,NA)
      Parent<-c(Parent,Scount)
      is.extant<-c(is.extant,TRUE)

      # get a new speciation waiting time, and include it in the vector
      WaitTimeS<-ifelse(pp(tNow)>0,
                        rexp_var(1,pp,tNow,tmax,pshape,
                                 ifelse(TS[Scount]<0,0,TS[Scount]), fast),Inf)
    }

    # reached the time of extinction
    tNow<-tExp

    # record extinction, and if species is extant make it more than tmax
    TE[Scount]<-ifelse(tNow<tmax,tNow,tmax+0.01)
    is.extant[Scount]<-ifelse(TE[Scount]>tmax,TRUE,FALSE)

    # next species
    Scount<-Scount+1
  }

  # now we invert TE and TS so time goes from tmax to 0
  TE <- tmax - TE
  TS <- tmax - TS

  return(list(TE=TE,TS=TS,PAR=Parent,EXTANT=is.extant))
}
