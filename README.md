# paleobuddy

`paleobuddy` is an R package to simulate species diversification and generate fossil records and phylogenies. While the literature on species birth-death simulators is extensive, including important software like [paleotree](https://github.com/dwbapst/paleotree) and [APE](https://github.com/cran/ape), we concluded there were interesting gaps to be filled regarding possible diversification scenarios. Differently from most simulators in the field, we strived for flexibility over focus, implementing and planning to implement a large array of regimens for users to experiment with and combine. In this way, `paleobuddy` can be used in complement to other simulators as a flexible jack of all trades, or, in the case of scenarios implemented only here, can allow for robust and easy simulations for novel situations.

The latest version can be installed in R using the `devtools` package

```
library(devtools)
install_github("brpetrucci/paleobuddy")
```

## Important functions

`bd.sim` is the main birth-death simulator of `PaleoBuddy`, allowing for multiple arguments to build a large number of possible scenarios. One can choose any type of time-varying constant for speciation `pp` and extinction `qq`. On top of the base rates, we allow for a `shape` parameter for each, if one chooses to interpret `pp` and `qq` as scales of a Weibull distribution for age-dependent diversification. We take the novel step of including time-varying scale for the Weibull distribution as an option for the rates. While time-varying shape is implemented, it has not been thoroughly tested as of the writing of this, and so we do not recommend its use unless the user tests it themselves. One can also supply an `env` parameter to make rates dependent on an environmental variable such as temperature. Finally, one could supply rates as a numeric vector, and supply a corresponding `shifts` with the respective shift times. These can all be combined as the user wishes, creating a myriad of possible scenarios we believe will allow for unprecedented flexibility in a researcher's simulation tools.

`sample.clade` is a fossil record-generating function, returning an organized data frame with occurrence times - or occurrence time ranges, provided the user supplies the respective interval vector. It allows for a sampling rate `rr` that can be similarly flexible to `pp` and `qq` above, with the exception of a `shape` parameter, since we ommitted that option given the absence of the use of Weibull distributions to model age-dependent sampling in the literature. Instead, we allow for the user to supply a function they wish to use as age-dependent sampling, `pFUN`, such as the PERT distribution used in [PyRate](https://github.com/dsilvestro/PyRate). If possible, the user can supply a maximizer for that function, `pFUNMax`, which would lead to faster computation. In the case of sampling, age-dependency and time-varying sampling rates are incompatible, at least as of the initial publication of the package. Still, `SampleClade` allows for unprecedented flexibility in sampling, letting the user combine as they wish time-varying, and environmentally-dependent functions, and any maximizable age-dependent function as a sampling rate.

`make.phylo` closes the trio of most important functions of the package, taking a PB simulation and returning a `phylo` object from the APE package (see above).

Besides its main species diversification-simulating functions, `PaleoBuddy` also supplies the user with a few interesting statistical tools, such as `rexp.var`, a generalization of the `rexp` function in BaseR that allows for time-varying exponential rates and a `shape` parameter, in which case it generalizes the `rweibull` function.

## Getting started

The user can simulate a group as follows

```
n0 <- 1 # initial number of species
pp <- 0.1 # speciation rate
qq <- 0.05 # extinction rate
tMax <- 10 # maximum simulation time
sim <- bd.sim(n0, pp, qq, tMax)
```

One can generate more complex simulations with time-varying or age-dependent rates as well

```
n0 <- 1 # initial number of species
pp <- function(t) 0.1 + 0.025*t # speciation rate
qq <- 10 # extinction rate scale
qShape <- 2 # extinction rate shape
tMax <- 10 # maximum simulation time
sim <- bd.sim(n0, pp, qq, tMax, qShape = qShape)
```

For details of further flexibility in rates, check `?bd.sim`.

The user can then use the simulated clade to generate fossil records or phylogenetic trees

```
rr <- 1 # sampling rate
bins <- seq(10, 0, -1) # something to simulate geologic intervals
samp <- sample.clade(sim = sim, rr = rr, tMax = tMax, bins = bins) # get a data frame with min and max ages of fossil occurrences
phy <- make.phylo(sim) # make a phylogenetic tree with the simulated group
ape::plot.phylo(phy) # plot it (requires APE)
ape::plot.phylo(ape::drop.fossil(phy)) # plot the molecular phylogeny
```

## Authors

`paleobuddy` was idealized by Bruno do Rosario Petrucci and Tiago Bosisio Quental. The birth-death, statistical and part of the sampling functions were written by Bruno. Most of the sampling functions were written by Matheus Januário.
