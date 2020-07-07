#' Returns the information for a simulated clade given constant speciation
#' and extinction rates
#'
#' \code{BDSimConstant} takes a rate, a maximum simulation time and an initial
#' number of species and returns the history of the clade originated from that
#' number of species with rates equal to the given rates. Time runs from
#' \code{tmax} to 0, to be consistent with the literature, though one can
#' easily invert that by subtracting the results from \code{tmax}.
#'
#' @param N0 initial number of species, usually 1. Good param to
#' tweak if one is observing a low sample size when testing.
#'
#' @param p speciation rate. Must be constant. No error check since
#' in the package this function will only be called by \code{BDSim}, which will
#' only call it if the rate is constant.
#'
#' @param q extinction rate, similar to above.
#'
#' @param tmax ending time of simulation. Any species still living
#' after \code{tmax} is considered extant, and any species that would be
#' generated after \code{tmax} is not born.
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
#' # first we define a function to calculate the mean diversity and var at time t
#' SimMean<-function(t, SimList){
#'   SimExtantT<-as.numeric(lapply(1:length(SimList),function(y){
#'     TS <- tmax - SimList[[y]]$TS
#'     TE <- tmax - SimList[[y]]$TE
#'     length(which(TS<=t&TE>=t))}))
#'   return(list(mean=mean(SimExtantT), var=(var(SimExtantT))))
#' }
#' # note the tmax -, rescaling the vectors so we can work only with functions
#' # going forward in time
#'
#' # also, we need functions to calculate the expected var at time t
#' Int<-function(t, div) {
#'   return(integrate(div, 0, t)$value)
#' }
#' DivVar<-function(t, div, qq){
#'   return(N0*exp(Int(t, div))*(exp(Int(t, div)) - 1 + 2*exp(Int(t,div))*
#'                                 integrate(Vectorize(function(x)
#'                                   exp(-Int(x, div))*qq(x)), 0, t)$value))
#' }
#' # and a time parameter we will need
#' tmax <- 40
#' Time <- 1:tmax
#'
#' # now we can test a couple scenarios
#' # first, extinction 0
#' N0 <- 1
#' tmax <- 40
#' p <- 0.1
#' q <- 0
#' SimList <- lapply(1:10000, function(x) BDSimConstant(N0, p, q, tmax))
#'
#' # let us make vectors to hold the average diversity and variance
#' div <- Vectorize(function(t) p - q)
#' qq <- Vectorize(function(t) q)
#'
#' MeanDiv <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$mean))
#' ExpectedDiv <- VarRateExp(div, 1, Time)
#'
#' MeanVar <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$var))
#' ExpectedVar <- unlist(lapply(Time, function(x) DivVar(x, div, qq)))
#'
#' # and now let us check out the plots
#' plot(Time, log(MeanDiv), type='l')
#' lines(Time, log(ExpectedDiv), col='RED')
#' plot(Time, log(MeanVar), type='l')
#' lines(Time, log(ExpectedVar), type='l', col='RED')
#'
#' # now let us try to turn extinction up a bit
#' N0 <- 1
#' tmax <- 40
#' p <- 0.1
#' q <- 0.04
#' SimList <- lapply(1:10000, function(x) BDSimConstant(N0, p, q, tmax))
#'
#' div <- Vectorize(function(t) p - q)
#' qq <- Vectorize(function(t) q)
#'
#' MeanDiv <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$mean))
#' ExpectedDiv <- VarRateExp(div, N0, Time)
#'
#' MeanVar <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$var))
#' ExpectedVar <- unlist(lapply(Time, function(x) DivVar(x, div, qq)))
#'
#' plot(Time, log(MeanDiv), type='l')
#' lines(Time, log(ExpectedDiv), col='RED')
#' plot(Time, log(MeanVar), type='l')
#' lines(Time, log(ExpectedVar), type='l', col='RED')
#'
#' # we can also try a pure-death process, starting with more species
#' N0 <- 100
#' tmax <- 40
#' p <- 0
#' q <- 0.02
#' SimList <- lapply(1:10000, function(x) BDSimConstant(N0, p, q, tmax))
#'
#' div <- Vectorize(function(t) p - q)
#' qq <- Vectorize(function(t) q)
#'
#' MeanDiv <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$mean))
#' ExpectedDiv <- VarRateExp(div, N0, Time)
#'
#' MeanVar <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$var))
#' ExpectedVar <- unlist(lapply(Time, function(x) DivVar(x, div, qq)))
#'
#' plot(Time, log(MeanDiv), type='l')
#' lines(Time, log(ExpectedDiv), col='RED')
#' plot(Time, log(MeanVar), type='l')
#' lines(Time, log(ExpectedVar), type='l', col='RED')
#'
#' # all the cases seem to agree pretty well with expectations
#'
#' @name BDSimConstant
#' @rdname BDSimConstant
#' @export


BDSimConstant<-function(N0 = 1, p, q, tmax){
  # initialize the vectors to hold times of speciation and extinction, parents
  # and status (extant or not)
  TS<-rep(-0.01,N0)
  TE<-rep(NA,N0)
  Parent<-rep(NA,N0)
  is.extant<-rep(TRUE,N0)

  # initialize the counting variable
  Scount<-1

  # while we have more species in a vector than we have analyzed,
  while (length(TE)>=Scount){
    # TS starts at -0.01 to show it was alive at the beginning, but to count
    # time we need to start at 0
    tNow<-ifelse(TS[Scount]<0,0,TS[Scount])

    # draw waiting times with rexp()
    WaitTimeS<-ifelse(p>0,rexp(1,p),Inf)
    WaitTimeE<-ifelse(q>0,rexp(1,q),Inf)

    # if the time of extinction is after the end of the simulation, make it tmax
    tExp<-min(tNow+WaitTimeE, tmax)

    # while there are fast enough speciations before the species goes extinct,
    while ((tNow+WaitTimeS)<=tExp){
      # update time
      tNow<-tNow+WaitTimeS

      # create a new species with corresponding TE, TS and parent
      TS<-c(TS,tNow)
      TE<-c(TE,NA)
      Parent<-c(Parent,Scount)
      is.extant<-c(is.extant,TRUE) # it is alive

      # take a new waiting time - if now + WaitTimeS is still less than when
      # the species goes extinct, repeat
      WaitTimeS<-ifelse(p>0,rexp(1,p),Inf)
    }

    # reached the time of the species extinction
    tNow<-tExp

    # record the extinction - if tExp >= tmax, it didn't go extinct
    TE[Scount]<-ifelse(tNow<tmax,tNow, tmax+0.01)
    is.extant[Scount]<-ifelse(TE[Scount] > tmax,TRUE,FALSE)

    # next species
    Scount<-Scount+1
  }

  # finally, we invert both TE and TS to attain to the convention that time
  # runs from 0 to tmax
  TE <- tmax - TE
  TS <- tmax - TS

  return(list(TE=TE,TS=TS,PAR=Parent,EXTANT=is.extant))
}
