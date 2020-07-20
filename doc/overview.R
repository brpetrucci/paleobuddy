## ----results = "hide"---------------------------------------------------------
# importing the package functions
library(paleobuddy)

## -----------------------------------------------------------------------------
# we set a seed so the results are reproducible
set.seed(1)

# set the necessary parameters
# initial number of species
n0 <- 1

# speciation rate - approx. 1 speciation event every 4my
# we are trying to create a big phylogeny so phytools can function better
pp <- 0.25

# extinction rate - approx. 1 extinction event every 10my
qq <- 0.15

# maximum simulation time - species that die after this are considered extant
tMax <- 50

# run the simulation
sim <- bd.sim(n0, pp, qq, tMax)

# take a look at the way the result is organized
sim

## -----------------------------------------------------------------------------
# there are currently not many customization options for phylogenies
phy <- make.phylo(sim)

# plot it with APE - hide tip labels since there are a lot so it looks cluttered
ape::plot.phylo(phy, show.tip.label = FALSE)
ape::axisPhylo()

# plot the molecular phylogeny
ape::plot.phylo(ape::drop.fossil(phy), show.tip.label = FALSE)
ape::axisPhylo()

## ----eval=FALSE---------------------------------------------------------------
#  # find the fit
#  fit <- phytools::fit.bd(ape::drop.fossil(phy))
#  
#  # get the birth and death rate
#  b <- fit$b
#  d <- fit$d

## -----------------------------------------------------------------------------
# set a seed
set.seed(2)

# create simulation
# note nFinal and extOnly, defining we want 200 ormore extant species at the end
sim <- bd.sim(n0, pp, qq, tMax, nFinal = c(200, Inf), extOnly = TRUE)

# check the number of extant species
sum(sim$EXTANT)

## ----eval=FALSE---------------------------------------------------------------
#  # find the fit
#  fit <- phytools::fit.bd(ape::drop.fossil(make.phylo(sim)))

## -----------------------------------------------------------------------------
# we set a seed so the results are reproducible
set.seed(5)

# set the necessary parameters
# initial number of species
n0 <- 1

# speciation rate - it can be any function of time!
pp <- function(t) {
  0.1 + 0.001*t
}

# extinction rate - also can be any function of time
qq <- function(t) {
  0.03 * exp(-0.01*t)
}

# maximum simulation time - species that die after this are considered extant
tMax <- 50

# run the simulation
sim <- bd.sim(n0, pp, qq, tMax)

# check the resulting clade out
ape::plot.phylo(make.phylo(sim), show.tip.label = FALSE)

## -----------------------------------------------------------------------------
# again set a seed
set.seed(1)

# set the sampling rate
# using a simple case - there will be on average T occurrences per species,
# where T is the species duration
rr <- 1

# run the sampling simulation only for the first 15 species for brevity's sake
samp <- suppressMessages(sample.clade(S = 1:15, sim = sim, rr = rr, tMax = tMax))
# suppressing messages - the message is to inform the user how many species
# left no fossils. In this case, it was 0

# take a look at how the output is organized
head(samp)

## -----------------------------------------------------------------------------
# make a copy
pSamp <- samp

# change the extant column
pSamp["Extant"][pSamp["Extant"] == FALSE] = "extant"
pSamp["Extant"][pSamp["Extant"] == TRUE] = "extinct"

# change column names
colnames(pSamp) <- c("Species", "Status", "min_age", "max_age")

# check it out
head(pSamp)

## -----------------------------------------------------------------------------
per.capita <-function(faBins, laBins, bins) {
  # create vectors to hold species that were born before and die in interval i,
  # species who were born in i and die later,
  # and species who were born before and die later
  NbL <- NFt <- Nbt <- rep(0, length(bins))
  
  # for each interval
  for (i in 1:length(bins)) {
    # number of species that were already around before i and are not seen again
    NbL[i] <- sum(faBins > bins[i] & laBins == bins[i])
    
    # number of species that were first seen in i and are seen later
    NFt[i] <- sum(faBins == bins[i] & laBins < bins[i])
    
    # number of species that were first seen before i and are seen after i
    Nbt[i] <- sum(faBins > bins[i] & laBins < bins[i])
  }
  
  # calculate the total rates
  p <- log((NFt + Nbt) / Nbt)
  q <- log((NbL + Nbt) / Nbt)
  return(list(p = p, q = q))
}

## -----------------------------------------------------------------------------
# get the species names
ids <- unique(samp$Species)

# get the first appearance bins - the first time in bins where the fossil was seen (lower bound)
faBins <- unlist(lapply(ids, function(i) max(samp$MaxT[samp$Species == i])))

# get the last appearance bins - last time in bins where the fossil was seen (upper bound)
laBins <- unlist(lapply(ids, function(i) min(samp$MinT[samp$Species == i])))

# create the bins vector we have been using
bins <- seq(tMax, 0, -0.1)
# note this has a high resolution, the actual stratigraphic ranges are much coarser

# get the estimates
pc <- per.capita(faBins, laBins, bins)

## -----------------------------------------------------------------------------
# set a seed
set.seed(4)

# parameters to set things up
n0 <- 1
tMax <- 20

# speciation can be dependent on an environmental variable as well as time
pp <- function(t, env) {
  0.01 * t + 0.1*exp(0.1*env)
}

# let us generate some fake environmental data - you can pretend this is 
# percentage of max temperature or something
Env <- data.frame(1:20, runif(20))
# the package RPANDA supplies some useful environmental data frames - see ?bd.sim and ?bd.sim general for
# examples with real data

# we can make extinction be age-dependent by creating a shape parameter
qq <- 10
qShape <- 2
# this will make it so species durations are distributed as a Weibull with scale 10 and shape 2

# run the simulation
# we pass the shape and environmental parameters
sim <- bd.sim(n0, pp, qq, tMax, qShape = qShape, envPP = Env)
# note that pShape and envQQ also exist
# the defaults for all of these customization options is NULL

# check out the phylogeny
ape::plot.phylo(make.phylo(sim), show.tip.label = FALSE)

## -----------------------------------------------------------------------------
# speciation may be a step function, presented as a rates vector and a shift times vector
pList <- c(0.1, 0.2, 0.05)
pShifts <- c(0, 10, 15)
# pShifts could also be c(tMax, tMax - 10, tMax - 15) for identical results

# make.rate is the function bd.sim calls to create a rate, but we will use it here to see how it looks
pp <- make.rate(pList, tMax = 20, fShifts = pShifts)
t <- seq(0, 20, 0.1)
plot(t, pp(t), type = 'l', main = "Step function speciation rate",
     xlab = "Time (My)", ylab = "Rate")

# it is not possible to combine the pp(t, env) and pList methods listed above
# but we can create a step function dependent on environmental data with ifelse
q <- function(t, env) {
  ifelse(t < 10, env,
         ifelse(t < 15, env * 2, env / 2))
}

# pass it to make.rate with environmental data
qq <- make.rate(q, envF = Env)
plot(t, qq(t), type = 'l', main = "Environmental step function extinction rate",
     xlab = "Time (My)", ylab = "Rate")

## ----eval=FALSE---------------------------------------------------------------
#  # set seed again
#  set.seed(1)
#  
#  # age-dependent speciation with a step function
#  pList <- c(10, 5, 12)
#  pShifts <- c(0, 10, 15)
#  pShape <- 0.5
#  
#  # age-dependent extinction with a step function of an environmental variable
#  q <- function(t, env) {
#    ifelse(t < 10, 20*env,
#           ifelse(t < 15, 10*env, 13*env / 2))
#  }
#  qShape <- 2
#  
#  # run the simulation
#  sim <- bd.sim(n0, pList, q, tMax, pShape = 0.5, qShape = 2, envQQ = Env, pShifts = pShifts)

## -----------------------------------------------------------------------------
# as an example, we will use a PERT distribution, a hat-shaped distribution used in PyRate
# preservation function
dPERT <- function(t, s, e, sp, a = 3, b = 3, log = FALSE) {

  # check if it is a valid PERT
  if (e >= s) {
    message("There is no PERT with e >= s")
    return(rep(NaN, times = length(t)))
  }

  # find the valid and invalid times
  id1 <- which(t <= e | t >= s)
  id2 <- which(!(t <= e | t >= s))
  t <- t[id2]

  # initialize result vector
  res <- vector()

  # if user wants a log function
  if (log) {
    # invalid times get -Inf
    res[id1] <- -Inf

    # valid times calculated with log
    res[id2] <- log(((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b)))
  }
  # otherwise
  else{
    res[id1] <- 0

    res[id2] <- ((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b))
  }

  return(res)
}

# function to calculate max of the PERT - makes the sampling faster
dPERTmax <- function(s, e, sp) {
  return(((s - e) / 2) + e)
}

# set seed
set.seed(1)

# generate a quick simulation
sim <- bd.sim(n0, pp = 0.1, qq = 0.05, tMax = 20)

# another seed
set.seed(1)

# sample - high sampling rate so we can visualize the age-dependency
samp <- sample.clade(sim = sim, rr = 30, tMax = 20, bins = bins,
                  dFun = dPERT, dFunMax = dPERTmax)

# extract species identity
ids <- unique(samp$Species)

# approximate sampling time (since it is a range)
mids <- (samp$MaxT - samp$MinT) / 2 + samp$MinT

# for each species
for (i in 1:length(ids)) {
  # get the species number
  sp <- unique(as.numeric(gsub("spp_", "", ids[i])))

  # check the histogram
  hist(mids[samp$Species == ids[[i]]],
       main = paste0("spp = ", sp, "; duration ~ ",
                     round(sim$TS[sp] - sim$TE[sp], digits = 2), "my"),
       xlab = "Time (My)", probability = TRUE)

  # expected curve
  curve(dPERT(x, s = sim$TS[sp], e = sim$TE[sp], sp = sp), from = sim$TE[sp],
        to = sim$TS[sp], add = TRUE, col = "red", n = 100)
}

## -----------------------------------------------------------------------------
# simple simulation, starting with more than one species
sim <- bd.sim(n0 = 2, pp = 0.1, qq = 0, tMax = 20, nFinal = c(20, Inf))

# separate the lineages
clades <- find.lineages(sim)

# plot each phylogeny

# clade 1
ape::plot.phylo(make.phylo(clades$clade_1), show.tip.label = FALSE)

# clade 2
ape::plot.phylo(make.phylo(clades$clade_2), show.tip.label = FALSE)

