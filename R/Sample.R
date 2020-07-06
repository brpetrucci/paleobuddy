#' Species sampling
#'
#' \code{Sample} takes a species number, a vector of speciation and extinction times,
#'  a sampling rate with possible shape and a maximum time for simulation.
#'
#' @param S the species number to be sampled. Since \code{Sample} will be
#' called by a wrapper using \code{lapply}, it is through \code{S}
#' that we apply this function.
#'
#' @param TE a vector of extinction times, usually an output of \code{BDSim}.
#'
#' @param TS a vector of speciation times, usually an output of \code{BDSim}.
#'
#' @param rr a sampling rate function. Can be created by \code{MakeRate} for
#' simplicity, but can be any time-varying function.
#'
#' @param tmax the maximum simulation time, used by \code{rexp_var}.
#'
#' @return a list of occurrences for that species.
#'
#' @author written by Bruno do Rosario Petrucci and Matheus Januario.
#'
#' @examples
#'
#' # Note: all examples use just 1 lineage and very large preservation rates just
#' # to be clearer to the reader. Most of times preservation will never be this high
#'
#' # let us start with a linear increase in preservation rate
#'
#' sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
#' while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived
#'                                  # lineage which will obscure the pattern
#'   sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
#' }
#'
#' # preservation function
#' r<-function(t) {
#'   return(50+t*30)
#' }
#'
#' t<-seq(0, 10, by=.1)
#' # visualizing from the past to the present
#' plot(x=t, y=rev(r(t)), main="Simulated preservation", type="l", col="red",
#' xlab="Mya", ylab="preservation rate", xlim=c(10, sim$TE[1]))
#'
#' occs<-Sample(S = 1, TS = sim$TS[1], TE = sim$TE[1], rr = r, tmax = 10)
#' hist(occs,
#'      xlim=c(10, sim$TE[1]), #changing axis
#'      xlab="Mya") #informative labels
#' lines(t, rev(r(t)))
#'
#' # now let us try a step function
#'
#' sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
#' while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived
#'                                  # lineage which will obscure the pattern
#'   sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
#' }
#'
#' # we can create the sampling rate here from a few vectors
#' rlist <- c(50, 20, 80)
#' rshifts <- c(0, 4, 8) # this could be c(10, 6, 2) and would produce the same function
#'
#' r <- MakeRate(rlist, tmax=10, fshifts=rshifts)
#'
#' t<-seq(0, 10, by=.1)
#' plot(x=t, y=rev(r(t)), main="Simulated preservation", type="l", col="red",
#' xlab="Mya", ylab="preservation rate", xlim=c(10, sim$TE[1]))
#'
#' occs<-Sample(S = 1, TS = sim$TS[1], TE = sim$TE[1], rr = r, tmax = 10)
#' hist(occs,
#'      xlim=c(10, sim$TE[1]), #changing axis
#'      xlab="Mya") #informative labels
#' abline(v=c(6,2), col="red") # frontiers of each regime
#'
#' # we can create a step function in a different way as well
#' sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
#' while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived lineage
#'                                  # which will obscure the pattern
#'   sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
#' }
#'
#' # preservation function
#' r<-function(t) {
#'   ifelse(t < 4, 50,# preservation rates at first regime (4 - 0 Mya)
#'          ifelse(t < 8, 20, # preservation rates at second regime (8 - 5 Mya)
#'                 80)) # preservation rates at third regime (after 8 Mya)
#' }
#'
#' t<-seq(0, 10, by=.1)
#' plot(x=t, y=rev(r(t)), main="Simulated preservation", type="l", col="red",
#' xlab="Mya", ylab="preservation rate", xlim=c(10, sim$TE[1]))
#'
#' occs<-Sample(S = 1, TS = sim$TS[1], TE = sim$TE[1], rr = r, tmax = 10)
#' hist(occs,
#'      xlim=c(10, sim$TE[1]), # changing axis
#'      xlab="Mya") # informative labels
#' abline(v=c(6, 2), col="red") # frontiers of each regime
#'
#' @name Sample
#' @rdname Sample
#' @export

Sample<-function(S,TE,TS,rr,tmax){
  TE <- tmax - TE
  TS <- tmax - TS

  # start when the species was born, end when it died
  Now<-ifelse(TS[S]<0,0,TS[S])
  End<-ifelse(TE[S]>tmax,tmax,TE[S])

  # make rr a function if it isn't
  r <- rr
  rr <- ifelse(is.numeric(r), Vectorize(function(t) r), r)

  # initialize vector
  sampled<-c()

  # while we have not reached the time of extinction
  while (Now<End){
    # take th waiting time for sampling, using rexp_var()
    WaitTimeR<-ifelse(rr(Now)>0,rexp_var(1,rr,Now,tmax), Inf)

    # advance in time
    Now<-Now+WaitTimeR

    # if sampling comes after extinction, we don't include this occurrence
    if (Now>End) break

    # add to the vector
    sampled<-c(sampled,Now)
  }

  # finally, we invert time so that it goes from tmax to 0
  sampled <- tmax - sampled

  return(sampled)
}
