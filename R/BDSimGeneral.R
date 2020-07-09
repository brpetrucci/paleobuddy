#' Non-constant rate Birth-Death simulation
#'
#' \code{BDSimGeneral} takes an initial number of species, speciation and
#' extinction rates (either constants, functions of time, or of time and an
#' environmental variable), a maximum simulation time and possibly a shape for
#' age-dependent speciation and/or extinction. It then generates the speciation 
#' and extinction times, parent and status information for the species in the 
#' simulation. Time runs from \code{tMax} to 0, to be consistent with the 
#' literature, though one can easily invert that by subtracting the results from
#' \code{tMax}.
#'
#' @param n0 initial number of species, usually 1. Good parameter
#' to tweak if one is observing a low sample size when testing.
#'
#' @param pp function to hold the speciation rate over time. It will either be
#' interpreted as an exponential rate, or a Weibull scale if 
#' \code{pShape != NULL}.
#'
#' @param qq similar to above, but for the extinction rate.
#' 
#' Note: this function is meant to be called by \code{BDSim}, so it neither
#' allows for as much flexibility, nor call \code{MakeRate}. If the user wishes
#' to use \code{BDSimGeneral} with environmental or step-function rates, they
#' can generate the rate with \code{MakeRate} and supply it to the function.
#'
#' @param tMax ending time of simulation. Any species still living
#' after \code{tMax} is considered extant, and any species that would be
#' generated after \code{tMax} is not born.
#'
#' @param pShape shape parameter for the Weibull distribution for age-dependent
#' speciation. Default is \code{NULL}, where \code{pp} will be considered a
#' time-dependent exponential rate. For \code{pShape != NULL}, \code{pp} will
#' be considered a scale, and \code{rexp_var} will draw a Weibull distribution
#' instead.
#'
#' @param qShape similar as above, but for the extinction rate.
#' 
#' @param fast used for \code{BDSimGeneral}. When \code{TRUE}, sets 
#' \code{rexp_var} to throw away waiting times higher than the maximum 
#' simulation time. Should be \code{FALSE} for unbiased testing of age 
#' dependency. User might also se it to \code{FALSE} for more accurate waiting
#' times.
#' 
#' @param trueExt used for \code{BDSimGeneral}. When \code{TRUE}, time of 
#' extinction of extant species will be the true time, otherwise it will be 
#' tMax+0.01. Need to be \code{TRUE} when testing age-dependent 
#' extinction.
#'
#' @return a list of vectors, as follows
#'
#' \describe{
#' \item{\code{TE}}{list of extinction times, with -0.01 as the time of
#' extinction for extant species.}
#'
#' \item{\code{TS}}{list of speciation times, with tMax+0.01 as the time of
#' speciation for species that started the simulation.}
#'
#' \item{\code{PAR}}{list of parents. Species that started the simulation have
#' NA, while species that were generated during the simulation have their
#' parent's number. Species are numbered as they are born.}
#'
#' \item{\code{EXTANT}}{list of booleans representing whether each species is
#' extant.}}
#'
#' @author written by Bruno do Rosario Petrucci.
#'
#' @examples
#'
#' # first we define a function to calculate the mean diversity and var at time t
#' sim.mean<-function(t, simList) {
#'   # apply to each replicate of the simulation
#'   SimExtantT<-as.numeric(lapply(1:length(simList), function(y) {
#'     
#'     # invert extinction and speciation times, so that we may consider functions
#'     # running from past to present
#'     TS <- tMax - simList[[y]]$TS
#'     TE <- tMax - simList[[y]]$TE
#'     
#'     # find how many species are alive at time t
#'     length(which(TS <= t & TE >= t))}))
#'   return(list(mean = mean(SimExtantT), var = (var(SimExtantT))))
#' }
#' 
#' # also, we need functions to calculate the expected varianc at time t
#' int<-function(t, div) {
#'   # integrate diversity from 0 to time t
#'   return(integrate(div, 0, t)$value)
#' }
#' 
#' # formula for the variance - see Kendall 1948
#' div.var<-function(t, div, qq) {
#'   return(n0*exp(int(t, div))*(exp(int(t, div)) - 1 + 2*exp(int(t, div))*
#'                                 integrate(Vectorize(function(x)
#'                                   exp(-int(x, div))*qq(x)), 0, t)$value))
#' }
#' 
#' # and a time parameter we will need for plotting
#' tMax <- 40
#' time <- 1:tMax
#' 
#' # now we can test a couple scenarios
#' # first, even though this is BDSimGeneral, we can try constant rates
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- 0.11
#' 
#' # extinction
#' q <- 0.08
#' 
#' \dontrun{
#' # run the simulations
#' simList <- lapply(1:10000, function(x) BDSimGeneral(n0, p, q, tMax))
#' 
#' # let us make vectors to hold the average diversity and variance
#' 
#' # function for diversity
#' div <- Vectorize(function(t) p - q)
#' 
#' # function for extinction - needed for the variance
#' qq <- Vectorize(function(t) q)
#' 
#' # calculate the mean diversity at our time points
#' MeanDiv <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$mean))
#' 
#' # calculate the expected diversity for the sime time points
#' expectedDiv <- VarRateExp(div, 1, time)
#' 
#' # do the same with variance
#' meanVar <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$var))
#' expectedVar <- unlist(lapply(time, function(x) div.var(x, div, qq)))
#' 
#' # and now let us check out the plots
#' 
#' # mean diversity, compared with expected
#' plot(time, log(MeanDiv), type = 'l', main = "Species diversity", 
#'      xlab = "Time (My)", ylab = "log(Diversity)")
#' lines(time, log(expectedDiv), col = 'RED')
#' legend(x = 5, y = log(max(MeanDiv)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' 
#' # same for variance
#' plot(time, log(meanVar), type = 'l',  main = "Species diversity variance", 
#'      xlab = "Time (My)", ylab = "log(Variance)")
#' lines(time, log(expectedVar), type = 'l', col='RED')
#' legend(x = 5, y = log(max(meanVar)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' }
#' 
#' # we can complicate things further with a linear function as a rate
#' # BDSimGeneral takes longer so we run examples for 1000 replicates instead
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' p <- function(t) {
#'   return(0.05 + 0.005*t)
#' }
#' 
#' # extinction
#' q <- 0.05
#' 
#' \dontrun{
#' # run the simulations
#' simList <- lapply(1:1000, function(x) BDSimGeneral(n0, p, q, tMax))
#' 
#' # let us make vectors to hold the average diversity and variance
#' 
#' # we won't need these for all simulations but it is good to make it uniform
#' # function for speciation
#' pp <- Vectorize(function(t) p(t))
#' 
#' # function for extinction
#' qq <- Vectorize(function(t) q)
#' 
#' # function for diversity
#' div <- Vectorize(function(t) p(t) - q)
#' 
#' # calculate the mean diversity at our time points
#' MeanDiv <- unlist(lapply(time, function(x) sim.mean(x, simList = simList)$mean))
#' 
#' # calculate the expected diversity for the sime time points
#' expectedDiv <- VarRateExp(div, 1, time)
#' 
#' # do the same with variance
#' meanVar <- unlist(lapply(time, function(x) sim.mean(x, SimList = SimList)$var))
#' expectedVar <- unlist(lapply(time, function(x) div.var(x, div, qq)))
#' 
#' # and now let us check out the plots
#' 
#' # mean diversity, compared with expected
#' plot(time, log(MeanDiv), type = 'l', main = "Species diversity", 
#'      xlab = "Time (My)", ylab = "log(Diversity)")
#' lines(time, log(expectedDiv), col = 'RED')
#' legend(x = 5, y = log(max(MeanDiv)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' 
#' # same for variance
#' plot(time, log(meanVar), type = 'l',  main = "Species diversity variance", 
#'      xlab = "Time (My)", ylab = "log(Variance)")
#' lines(time, log(expectedVar), type = 'l', col='RED')
#' legend(x = 5, y = log(max(meanVar)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' }
#' 
#' # we can also create a step function. Keep in mind this is a slower way than by
#' # creating step functions using ifelse()
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation rate
#' p <- function(t) {
#'   return(0.05 + 0.005*t)
#' }
#' 
#' # list of extinction rates
#' qList <- c(0.04, 0.06, 0.07)
#' 
#' # list of shift times. Note qShifts could be c(40, 20, 10) for 
#' # identical results
#' qShifts <- c(0, 20, 30)
#' 
#' # let us take a look at how MakeRate will make it a step function
#' q <- MakeRate(qList, fShifts = qShifts)
#' 
#' # and plot it
#' plot(seq(0, tMax, 0.1), q(seq(0, tMax, 0.1)), type = 'l',
#'      main = "Extintion rate as a step function", xlab = "Time (My)",
#'      ylab = "Rate (species/My)")
#' # note that this is slower than creating a step function with ifelse(), in this
#' # case q <- function(t) ifelse(t < 20, 0.04, ifelse(t < 30, 0.06, 0.07))
#' 
#' # also note that if done with ifelse(), the function must go from 0, instead of
#' # from tMax
#' 
#' \dontrun{
#' # run the simulations
#' SimList <- lapply(1:1000, function(x) BDSimGeneral(n0, p, q, tMax))
#' 
#' # let us make vectors to hold the average diversity and variance
#' 
#' # we won't need these for all simulations but it is good to make it uniform
#' # function for speciation
#' pp <- Vectorize(function(t) p(t))
#' 
#' # function for extinction
#' qq <- Vectorize(function(t) q(t))
#' 
#' # function for diversity
#' div <- Vectorize(function(t) p(t) - q(t))
#' 
#' # calculate the mean diversity at our time points
#' MeanDiv <- unlist(lapply(time, function(x) sim.mean(x, SimList = SimList)$mean))
#' 
#' # calculate the expected diversity for the sime time points
#' expectedDiv <- VarRateExp(div, 1, time)
#' 
#' # do the same with variance
#' meanVar <- unlist(lapply(time, function(x) sim.mean(x, SimList = SimList)$var))
#' expectedVar <- unlist(lapply(time, function(x) div.var(x, div, qq)))
#' 
#' # and now let us check out the plots
#' 
#' # mean diversity, compared with expected
#' plot(time, log(MeanDiv), type = 'l', main = "Species diversity", 
#'      xlab = "Time (My)", ylab = "log(Diversity)")
#' lines(time, log(expectedDiv), col = 'RED')
#' legend(x = 5, y = log(max(MeanDiv)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' 
#' # same for variance
#' plot(time, log(meanVar), type = 'l',  main = "Species diversity variance", 
#'      xlab = "Time (My)", ylab = "log(Variance)")
#' lines(time, log(expectedVar), type = 'l', col='RED')
#' legend(x = 5, y = log(max(meanVar)), legend = c("Expected", "Observed"),
#'        col = c("RED", "BLACK"), lty = c(1,1))
#' }
#' 
#' # another feature to add is age dependency. Note that since there is no 
#' # analytical solution to this system, we must test species longevity directly 
#' # and therefore must pass fast = FALSE to the function
#' 
#' if (requireNamespace("fitdistrplus", quietly = TRUE)) {
#'   # initial number of species
#'   n0 <- 1
#'   
#'   # maximum simulation time
#'   tMax <- 40
#'   
#'   # speciation
#'   p <- 0.15
#'   
#'   # extinction - a Weibull scale
#'   q <- 10
#'   
#'   # extinction shape
#'   qShape <- 1
#'   
#'   \dontrun{
#'   # run simulations - note fast = FALSE and trueExt = TRUE so we can accurately
#'   # fit the results to a Weibull
#'   SimList <- lapply(1:1000, function(x)
#'     BDSimGeneral(n0, p, q, tMax, qShape = qShape, fast = FALSE, trueExt = TRUE))
#'   
#'   # now we can use fitdistrplus to check that, on average, the longevities 
#'   # simulated follow a Weibull distribution
#'   
#'   # vectors to hold the shape and scale values
#'   shapes <- c()
#'   scales <- c()
#'   
#'   # for each simulation
#'   for (i in 1:length(SimList)) {
#'     
#'     # invert the time so we can perform a fit
#'     TE <- tMax - SimList[[i]]$TE
#'     TS <- tMax - SimList[[i]]$TS
#'     
#'     # TS could be -0.01, set it to 0 if so
#'     TS <- ifelse(TS < 0, 0, TS)
#'     
#'     # cannot fit for simulations that had only one species
#'     if (length(TE) < 2) next
#'     
#'     estimate <- fitdistrplus::fitdist(TE - TS, distr = "weibull", method = "mge",
#'                                       lower = c(0,0), start = list(shape = 1, scale = 10),
#'                                       gof = 'CvM')$estimate
#'     shapes <- c(shapes, estimate[1])
#'     scales <- c(scales, estimate[2])
#'   }
#'   
#'   # make a boxplot
#'   boxplot(shapes, outline = FALSE, main = "Boxplot of shapes")
#'   abline(h = qShape)
#'   boxplot(scales, outline = FALSE, main = "Boxplot of scales")
#'   abline(h = q)
#'   }
#' }
#' 
#' # we can also have time-varying scale, but the complexity of the system makes
#' # the only test possible be a direct calculation of the expected longevity.
#' # Since we have done that in the tests for rexp_var(), we will not repeat it
#' # here.
#' 
#' \dontrun{
#' # finally, we could have environmental dependency on a rate. To test that, we 
#' # need RPANDA
#' if (requireNamespace("RPANDA", quietly = TRUE)) {
#'   # initial number of species
#'   n0 <- 1
#'   
#'   # maximum simulation time
#'   tMax <- 40
#'   
#'   # temperature-dependent speciation
#'   p_t <- function(t, temp) {
#'     return(0.04*exp(0.15*temp))
#'   }
#'   
#'   # extinction
#'   q <- 0.01
#'   
#'   # using RPANDA to get the temperature data
#'   data(InfTemp, package="RPANDA")
#'   
#'   # speciation
#'   p <- MakeRate(p_t, tMax, envF = InfTemp)
#'   
#'   # since we need many species to be able to test this effectively using
#'   # RPANDA, and the rates become really noisy with temperature, we set
#'   # only 100 simulations to finish it in a reasonable time
#'   
#'   # run simulations
#'   SimList <- lapply(1:1000, function(x) BDSimGeneral(n0, p, q, tMax))
#'   
#'   # let us make vectors to hold the average diversity and variance
#'   
#'   # we won't need these for all simulations but it is good to make it uniform
#'   # function for speciation
#'   pp <- Vectorize(function(t) p(t))
#'   
#'   # function for extinction
#'   qq <- Vectorize(function(t) q)
#'   
#'   # function for diversity
#'   div <- Vectorize(function(t) p(t) - q)
#'   
#'   # calculate the mean diversity at our time points
#'   MeanDiv <- unlist(lapply(time, function(x) sim.mean(x, SimList = SimList)$mean))
#'   
#'   # calculate the expected diversity for the sime time points
#'   expectedDiv <- VarRateExp(div, 1, time)
#'   
#'   # do the same with variance
#'   meanVar <- unlist(lapply(time, function(x) sim.mean(x, SimList = SimList)$var))
#'   expectedVar <- unlist(lapply(time, function(x) div.var(x, div, qq)))
#'   
#'   # and now let us check out the plots
#'   
#'   # mean diversity, compared with expected
#'   plot(time, log(MeanDiv), type = 'l', main = "Species diversity", 
#'        xlab = "Time (My)", ylab = "log(Diversity)")
#'   lines(time, log(expectedDiv), col = 'RED')
#'   legend(x = 5, y = log(max(MeanDiv)), legend = c("Expected", "Observed"),
#'          col = c("RED", "BLACK"), lty = c(1,1))
#'   
#'   # same for variance
#'   plot(Time, log(meanVar), type = 'l',  main = "Species diversity variance", 
#'        xlab = "Time (My)", ylab = "log(Variance)")
#'   lines(Time, log(expectedVar), type = 'l', col='RED')
#'   legend(x = 5, y = log(max(meanVar)), legend = c("Expected", "Observed"),
#'          col = c("RED", "BLACK"), lty = c(1,1))
#'   
#'   # the plots look good, but the noise makes it confusing. We can also use the
#'   # fit_env function from RPANDA to test
#'   
#'   # create the necessary parameters for fit_env
#'   
#'   # vectors to hold the parameters estimated by fit_env - the value multiplying
#'   # the exponential and the value multiplying temperature
#'   intercepts <- c()
#'   mults <- c()
#'   
#'   # function fit_env will use to fit speciation
#'   f.l <- function(t, x, y) {
#'     y[1] * exp(y[2] * x)
#'   }
#'   
#'   # starting values for the fit
#'   lpar <- c(0.05, 0.2)
#'   
#'   # we fix mu since fit_env works much better this way. We could also set f.m 
#'   # to y[1] and set cst.mu = TRUE in fit_env, or not set any other options, but 
#'   # this leads to a much worse fit on the part of RPANDA
#'   f.m <- function(t, x, y) {
#'     0.06
#'   }
#'   
#'   # no starting values since there is nothing to estimate
#'   mpar <- c()
#'   
#'   # degrees of freedom for the data set
#'   dof <- smooth.spline(InfTemp[,1], InfTemp[,2])$df
#'   
#'   # matrix to hold estimated parameters
#'   par_matrix <- matrix(0, nrow = 0, ncol = 2)
#'   
#'   # for each replicate
#'   for (i in 1:length(SimList)) {
#'     # get the simulation in question
#'     sim <- SimList[[i]]
#'     
#'     # need more than 2 species to use fit_env
#'     if (length(sim$TE[sim$TE < 0]) < 2) next
#'     
#'     # get the phylogeny (needs to be molecular, so we use drop.fossil)
#'     phy <- ape::drop.fossil(MakePhylo(sim))
#'     
#'     # total time of the phylogeny
#'     # note picante is in RPANDA's namespace
#'     tot_time <- max(picante::node.age(phy)$ages)
#'     
#'     # perform the fit
#'     envFit <- RPANDA::fit_env(phy, InfTemp, tot_time, f.l, f.m, lpar,
#'                               mpar, df = dof, dt = 1e-3, fix.mu = TRUE)
#'     
#'     # append the estimated parameters to the matrix
#'     par_matrix <- rbind(par_matrix, envFit$lamb_par)
#'   }
#'   
#'   # one must remember that functions in RPANDA go from present to past, as
#'   # opposed to our functions. So we can test by comparing the value of our
#'   # rates with the estimates at the given time points
#'   
#'   # first get temperature temp at the time points
#'   find_t <- function(A, t) {
#'     # annoying expression but all it does is find the time point closest to t
#'     # in a vector A
#'     return(min(which(A == A[A - t == min(A[which(A - t > 0)] - t)])))
#'   }
#'   
#'   # find the corresponding Time points in InfTemp
#'   temp_t <- unlist(lapply(Time, function(x) find_t(InfTemp[, 1], x)))
#'   
#'   # get the temperature at those time points
#'   temp <- InfTemp[temp_t, 2]
#'   
#'   # get the value of the estimate and observed functions at these points
#'   estimate <- rev(mean(par_matrix[, 1]) *
#'                     exp(mean(par_matrix[, 2]) * temp))
#'   actual <- p_t(temp_t, temp)
#'   
#'   # total quadratic error should be low
#'   sum((estimate-actual)^2)
#' }
#' }
#' @name BDSimGeneral
#' @rdname BDSimGeneral
#' @export
#' 

BDSimGeneral<-function(n0,pp,qq,tMax,pShape=NULL,qShape=NULL,fast=TRUE,trueExt=FALSE) {
  # create vectors to hold times of speciation, extinction, parents and status
  TS<-rep(-0.01,n0)
  TE<-rep(NA,n0)
  parent<-rep(NA,n0)
  isExtant<-rep(TRUE,n0)

  # initialize species count
  sCount<-1

  # if shape is not null, make scale a function if it is not
  if (!is.null(pShape)) {
    p <- pp
    pp <- ifelse(is.numeric(p), Vectorize(function(t) p),
                 p)
  }
  if (!is.null(qShape)) {
    q <- qq
    qq <- ifelse(is.numeric(q), Vectorize(function(t) q),
                 q)
  }

  # while we have species to be analyzed still
  while (length(TE)>=sCount) {

    # get the time of speciation, or 0 if the species
    # was there at the beginning
    tNow<-ifelse(TS[sCount]<0,0,TS[sCount])

    # find the waiting time using rexp_var - note that in rexp_var we only
    # count t from tNow (to consider the rates as functions), so that
    # now we need to subtract tNow
    waitTimeS<-ifelse(is.numeric(pp), rexp(1, pp),
                      ifelse(pp(tNow)>0,
                             rexp_var(1,pp,tNow,tMax,pShape,
                                      ifelse(TS[sCount]<0,0,TS[sCount]), fast),Inf))
    waitTimeE<-ifelse(is.numeric(qq), rexp(1, qq),
                      ifelse(qq(tNow)>0,
                             rexp_var(1,qq,tNow,tMax,qShape,
                                      ifelse(TS[sCount]<0,0,TS[sCount]), fast),Inf))

    tExp<-tNow+waitTimeE

    # while there are fast enough speciations before the species goes extinct,
    while ((tNow+waitTimeS)<=min(tExp, tMax)) {

      # advance to the time of speciation
      tNow<-tNow+waitTimeS

      # add new times to the vectors
      TS<-c(TS,tNow)
      TE<-c(TE,NA)
      parent<-c(parent,sCount)
      isExtant<-c(isExtant,TRUE)

      # get a new speciation waiting time, and include it in the vector
      waitTimeS<-ifelse(is.numeric(pp), rexp(1, pp),
                        ifelse(pp(tNow)>0,
                               rexp_var(1,pp,tNow,tMax,pShape,
                                        ifelse(TS[sCount]<0,0,TS[sCount]), fast),Inf))
    }

    # reached the time of extinction
    tNow<-tExp

    # record extinction, and if species is extant make it more than tMax
    TE[sCount]<-ifelse(tNow<tMax | trueExt,tNow,tMax+0.01)
    isExtant[sCount]<-ifelse(TE[sCount]>tMax,TRUE,FALSE)

    # next species
    sCount<-sCount+1
  }

  # now we invert TE and TS so time goes from tMax to 0
  TE <- tMax - TE
  TS <- tMax - TS

  return(list(TE=TE,TS=TS,PAR=parent,EXTANT=isExtant))
}
