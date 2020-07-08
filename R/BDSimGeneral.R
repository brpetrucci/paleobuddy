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
#' testing of age dependency. User might also se it to \code{FALSE} for more
#' accurate waiting times.
#'
#' @param trueExt when \code{TRUE}, time of extinction of extant species will be
#' the true time, otherwise it will be tmax+0.01. Need to be \code{TRUE} when
#' testing age-dependent extinction
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
#'
#' # first, even though this is BDSimGeneral, we can try constant rates
#' N0 <- 1
#' tmax <- 40
#' p <- 0.11
#' q <- 0.08
#' SimList <- lapply(1:10000, function(x) BDSimGeneral(N0, p, q, tmax))
#'
#' # let us make vectors to hold the average diversity and variance
#' pp <- Vectorize(function(t) p)
#' qq <- Vectorize(function(t) q)
#' div <- Vectorize(function(t) pp(t) - qq(t))
#'
#' MeanDiv <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$mean))
#' ExpectedDiv <- VarRateExp(div, 1, Time)
#'
#' MeanVar <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$var))
#' ExpectedVar <- unlist(lapply(Time, function(x) DivVar(x, div, qq)))
#'
#' # and now let us check out the plots
#' plot(Time, log(MeanDiv), type='l', main="Species Diversity", xlab="Time (My)",
#'      ylab="log(Diversity)")
#' lines(Time, log(ExpectedDiv), col='RED')
#' legend(x=5, y=log(max(MeanDiv)), legend=c("Expected", "Observed"),
#'        col=c("RED", "BLACK"), lty=c(1,1))
#' plot(Time, log(MeanVar), type='l')
#' lines(Time, log(ExpectedVar), type='l', col='RED')
#' legend(x=5, y=log(max(MeanVar)), legend=c("Expected", "Observed"),
#'        col=c("RED", "BLACK"), lty=c(1,1))
#'
#' # we can complicate things a bit by making speciation time dependent
#' # note we lower the number of simulations, since BDSimGeneral takes
#' # longer for non-constant rates. This takes approximately 3 minutes,
#' # so if a user has time they can replicate this for more replicates
#' # and see it agrees with expectation even better
#' N0 <- 1
#' tmax <- 40
#' p <- function(t) {
#'   return(0.05+0.005*t)
#' }
#' q <- 0.05
#' SimList <- lapply(1:1000, function(x) BDSimGeneral(N0, p, q, tmax))
#'
#' # let us make vectors to hold the average diversity and variance
#' pp <- Vectorize(function(t) p(t))
#' qq <- Vectorize(function(t) q)
#' div <- Vectorize(function(t) pp(t) - qq(t))
#'
#' MeanDiv <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$mean))
#' ExpectedDiv <- VarRateExp(div, 1, Time)
#'
#' MeanVar <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$var))
#' ExpectedVar <- unlist(lapply(Time, function(x) DivVar(x, div, qq)))
#'
#' # and now let us check out the plots
#' plot(Time, log(MeanDiv), type='l', main="Species Diversity", xlab="Time (My)",
#'      ylab="log(Diversity)")
#' lines(Time, log(ExpectedDiv), col='RED')
#' legend(x=5, y=log(max(MeanDiv)), legend=c("Expected", "Observed"),
#'        col=c("RED", "BLACK"), lty=c(1,1))
#' plot(Time, log(MeanVar), type='l')
#' lines(Time, log(ExpectedVar), type='l', col='RED')
#' legend(x=5, y=log(max(MeanVar)), legend=c("Expected", "Observed"),
#'        col=c("RED", "BLACK"), lty=c(1,1))
#'
#' # we can also create a step function. Keep in mind this is a slower way than by
#' # creating step functions using ifelse()
#' N0 <- 1
#' tmax <- 40
#' p <- function(t) {
#'   return(0.05+0.005*t)
#' }
#' qlist <- c(0.04, 0.06, 0.07)
#' qshifts <- c(0, 20, 30)
#' q <- MakeRate(qlist, tmax, fshifts=qshifts)
#' SimList <- lapply(1:1000, function(x) BDSimGeneral(N0, p, q, tmax))
#'
#' # let us make vectors to hold the average diversity and variance
#' pp <- Vectorize(function(t) p(t))
#' qq <- Vectorize(function(t) q(t))
#' div <- Vectorize(function(t) pp(t) - qq(t))
#'
#' MeanDiv <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$mean))
#' ExpectedDiv <- VarRateExp(div, 1, Time)
#'
#' MeanVar <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$var))
#' ExpectedVar <- unlist(lapply(Time, function(x) DivVar(x, div, qq)))
#'
#' # and now let us check out the plots
#' plot(Time, log(MeanDiv), type='l', main="Species Diversity", xlab="Time (My)",
#'      ylab="log(Diversity)")
#' lines(Time, log(ExpectedDiv), col='RED')
#' legend(x=5, y=log(max(MeanDiv)), legend=c("Expected", "Observed"),
#'        col=c("RED", "BLACK"), lty=c(1,1))
#' plot(Time, log(MeanVar), type='l')
#' lines(Time, log(ExpectedVar), type='l', col='RED')
#' legend(x=5, y=log(max(MeanVar)), legend=c("Expected", "Observed"),
#'        col=c("RED", "BLACK"), lty=c(1,1))
#'
#' # another feature to add is age dependency. Note that since there is no analytical
#' # solution to this system, we must test species longevity directly and therefore
#' # must pass fast=FALSE to the function
#' N0 <- 1
#' tmax <- 40
#' p <- 0.15
#' q <- 10
#' qshape <- 1
#' SimList <- lapply(1:1000, function(x)
#'   BDSimGeneral(N0, p, q, tmax, qshape=qshape, fast=FALSE, trueExt=TRUE))
#'
#' # now we can use fitdistrplus to check that, on average, the longevities simulated
#' # follow a Weibull distribution
#' shapes <- c()
#' scales <- c()
#' for (i in 1:length(SimList)) {
#'   TE <- tmax - SimList[[i]]$TE
#'   TS <- tmax - SimList[[i]]$TS
#'   TS <- ifelse(TS<0, 0, TS)
#'
#'   if (length(TE) < 2) next
#'
#'   estimate <- fitdistrplus::fitdist(TE-TS,distr="weibull",method="mge",
#'                                     lower=c(0,0), start=list(shape=1,scale=10),
#'                                     gof='CvM')$estimate
#'   shapes <- c(shapes, estimate[1])
#'   scales <- c(scales, estimate[2])
#' }
#'
#' # make a boxplot
#' boxplot(shapes, outline=FALSE, main="Boxplot of shapes")
#' abline(h=qshape)
#' boxplot(scales, outline=FALSE, main="Boxplot of scales")
#' abline(h=q)
#'
#' # we can also have time-varying scale, but the complexity of the system makes
#' # the only test possible be a direct calculation of the expected longevity.
#' # Since we have done that in the tests for rexp_var(), we will not repeat it
#' # here.
#'
#' # finally, we could have environmental dependency on a rate. For that, we need
#' # RPANDA
#' if (requireNamespace("RPANDA", quietly=TRUE) &
#'     requireNamespace("ape", quietly=TRUE)) {
#'   N0 <- 1
#'   tmax <- 40
#'   p_t <- function(t, temp) {
#'     return(0.04*exp(0.15*temp))
#'   }
#'   q <- 0.01
#'
#'   # using RPANDA to get the temperature data
#'   data(InfTemp, package="RPANDA")
#'
#'   p <- MakeRate(p_t, tmax, env_f=InfTemp)
#'   # since we need many species to be able to test this effectively using
#'   # RPANDA, and the rates become really noisy with temperature, we set
#'   # only 100 simulations to finish it in a reasonable time
#'   SimList <- lapply(1:100, function(x) BDSimGeneral(N0, p, q, tmax))
#'
#'   # let us make vectors to hold the average diversity and variance
#'   pp <- Vectorize(function(t) p(t))
#'   qq <- Vectorize(function(t) q)
#'   div <- Vectorize(function(t) pp(t) - qq(t))
#'
#'   MeanDiv <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$mean))
#'   ExpectedDiv <- VarRateExp(div, 1, Time)
#'
#'   MeanVar <- unlist(lapply(Time, function(x) SimMean(x, SimList=SimList)$var))
#'   ExpectedVar <- unlist(lapply(Time, function(x) DivVar(x, div, qq)))
#'
#'   # and now let us check out the plots
#'   plot(Time, log(MeanDiv), type='l', main="Species Diversity", xlab="Time (My)",
#'        ylab="log(Diversity)")
#'   lines(Time, log(ExpectedDiv), col='RED')
#'   legend(x=5, y=log(max(MeanDiv)), legend=c("Expected", "Observed"),
#'          col=c("RED", "BLACK"), lty=c(1,1))
#'   plot(Time, log(MeanVar), type='l')
#'   lines(Time, log(ExpectedVar), type='l', col='RED')
#'   legend(x=5, y=log(max(MeanVar)), legend=c("Expected", "Observed"),
#'          col=c("RED", "BLACK"), lty=c(1,1))
#'
#'   # the plots look good, but the noise makes it confusing. We can also use the
#'   # fit_env function from RPANDA to test
#'   intercepts <- c()
#'   mults <- c()
#'
#'   # create the necessary parameters for fit_env
#'   f.l <- function(t, x, y) {
#'     y[1] * exp(y[2] * x)
#'   }
#'   lpar <- c(0.05, 0.2)
#'   # we fix mu since fit_env works much better this way. We could also set f.m to y[1]
#'   # and set cst.mu=TRUE in fit_env, or not set any other options, but this leads to
#'   # a much worse fit on the part of RPANDA
#'   f.m <- function(t, x, y) {
#'     0.06
#'   }
#'   mpar <- c()
#'
#'   dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df
#'
#'   par_matrix <- matrix(0, nrow = 0, ncol = 2)
#'   for (i in 1:length(SimList)) {
#'     sim <- SimList[[i]]
#'     if (length(sim$TE[sim$TE < 0]) < 2) next
#'
#'     # needs to be a molecular phylogeny
#'     phy <- ape::drop.fossil(MakePhylo(sim))
#'
#'     tot_time <- max(ape::node.age(phy)$ages)
#'     env_fit <- RPANDA::fit_env(phy,InfTemp,tot_time,f.l,f.m,lpar,mpar,df=dof,dt=1e-3, fix.mu=TRUE)
#'     par_matrix <- rbind(par_matrix, env_fit$lamb_par)
#'   }
#'
#'   # one must remember that functions in RPANDA go from present to past, as opposed
#'   # to our functions. So we can test by comparing the value of our rates with the
#'   # estimates at the given time points
#'
#'   # first get tempearature temp at the time points
#'   find_t <- function(A, t) {
#'    return(min(which(A==A[A-t==min(A[which(A-t>0)]-t)])))
#'   }
#'   temp_t <- unlist(lapply(Time, function(x) find_t(InfTemp[,1], x)))
#'   temp <- InfTemp[temp_t, 2]
#'
#'   # then get the value of the estimate and observed functions at these points
#'   estimate <- rev(mean(par_matrix[,1]) *
#'                     exp(mean(par_matrix[,2]) * temp))
#'   actual <- p_t(temp_t, temp)
#'
#'   # total quadratic error should be low
#'   sum((estimate-actual)^2)
#' }
#' @name BDSimGeneral
#' @rdname BDSimGeneral
#' @export

BDSimGeneral<-function(N0,pp,qq,tmax,pshape=NULL,qshape=NULL,fast=TRUE,trueExt=FALSE){
  # create vectors to hold times of speciation, extinction, parents and status
  TS<-rep(-0.01,N0)
  TE<-rep(NA,N0)
  Parent<-rep(NA,N0)
  is.extant<-rep(TRUE,N0)

  # initialize species count
  Scount<-1

  # if shape is not null, make scale a function if it is not
  if (!is.null(pshape)) {
    p <- pp
    pp <- ifelse(is.numeric(p), Vectorize(function(t) p),
                 p)
  }
  if (!is.null(qshape)) {
    q <- qq
    qq <- ifelse(is.numeric(q), Vectorize(function(t) q),
                 q)
  }

  # while we have species to be analyzed still
  while (length(TE)>=Scount){

    # get the time of speciation, or 0 if the species
    # was there at the beginning
    tNow<-ifelse(TS[Scount]<0,0,TS[Scount])

    # find the waiting time using rexp_var - note that in rexp_var we only
    # count t from tNow (to consider the rates as functions), so that
    # now we need to subtract tNow
    WaitTimeS<-ifelse(is.numeric(pp), rexp(1, pp),
                      ifelse(pp(tNow)>0,
                             rexp_var(1,pp,tNow,tmax,pshape,
                                      ifelse(TS[Scount]<0,0,TS[Scount]), fast),Inf))
    WaitTimeE<-ifelse(is.numeric(qq), rexp(1, qq),
                      ifelse(qq(tNow)>0,
                             rexp_var(1,qq,tNow,tmax,qshape,
                                      ifelse(TS[Scount]<0,0,TS[Scount]), fast),Inf))

    tExp<-tNow+WaitTimeE

    # while there are fast enough speciations before the species goes extinct,
    while ((tNow+WaitTimeS)<=min(tExp, tmax)){

      # advance to the time of speciation
      tNow<-tNow+WaitTimeS

      # add new times to the vectors
      TS<-c(TS,tNow)
      TE<-c(TE,NA)
      Parent<-c(Parent,Scount)
      is.extant<-c(is.extant,TRUE)

      # get a new speciation waiting time, and include it in the vector
      WaitTimeS<-ifelse(is.numeric(pp), rexp(1, pp),
                        ifelse(pp(tNow)>0,
                               rexp_var(1,pp,tNow,tmax,pshape,
                                        ifelse(TS[Scount]<0,0,TS[Scount]), fast),Inf))
    }

    # reached the time of extinction
    tNow<-tExp

    # record extinction, and if species is extant make it more than tmax
    TE[Scount]<-ifelse(tNow<tmax | trueExt,tNow,tmax+0.01)
    is.extant[Scount]<-ifelse(TE[Scount]>tmax,TRUE,FALSE)

    # next species
    Scount<-Scount+1
  }

  # now we invert TE and TS so time goes from tmax to 0
  TE <- tmax - TE
  TS <- tmax - TS

  return(list(TE=TE,TS=TS,PAR=Parent,EXTANT=is.extant))
}
