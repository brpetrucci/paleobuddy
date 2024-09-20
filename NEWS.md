<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# paleobuddy 1.1.0.0

## Simulating to a number of extant species

- `bd.sim` now allows for simulation up to a number of species as well
  as the original implementation simulating up to a certain time. To
  choose which realization of a certain number of species we choose, we
  use the method proposed by Stadler 2011 (see `bd.sim` documentation
  for full reference), whereby we simulate to a much larger number of
  species, then go back to check which periods had the desired number,
  and choose one through weighted sampling using the amount of time
  spent in each as the weights. Finally, we sample a uniform
  distribution to decide the length of the final period between when the
  species number was achieve and the end of the simulation.

## Trait-dependent dynamics

Added functions to simulate trait evolution and trait-dependent birth
death models, in particular State Speciation and Extinction (SSE) models
BiSSE, MuSSE, and QuaSSE. See documentation for full reference.

- `bd.sim.traits` simulates species diversification following a MuSSE
  model, where traits evolve from an Mk model and change speciation
  and/or extinction in a discrete fashion. Allows for the simulation of
  multiple traits, though currently the rates can only depend on one of
  them (each rate can depend on a different trait, though). Can simulate
  HiSSE as well, by setting the `nHidden` parameter, representing the
  number of hidden states, to something higher than 1. Can set separate
  number of states (observed or hidden), initial trait values, and
  transition matrices for each trait. See `?bd.sim.traits` for details.
  In the future, other trait-dependent diversification models will be
  implemented, including relaxing the assumption of discrete traits.
- `sample.clade.traits` simulates state-dependent fossil sampling,
  following a similar algorithm as `bd.sim.traits`. It allows for a
  hidden state model as well, but note that trait information (usually
  coming from `bd.sim.traits`) needs to include all states as observed
  for that (see examples in `?sample.clade.traits` and vignettes for
  details).
- `draw.sim` can now color longevity segments and fossil occurrences
  based on the trait values of a given trait for each species through
  time. Customization options include which colors to use for each
  state, whether to plot fossils as true occurrence times or ranges
  (already present in 1.0, but now with a dedicated argument for such,
  `fossilsFormat`), and where to place trait value legend. See
  `?draw.sim`, `?bd.sim.traits`, `?sample.clade.traits`, and the
  `overview` vignette for examples.

## Adding sampled ancestors to a phylogenetic tree

- `make.phylo` now allows for an optional `fossils` input representing a
  fossil record, and will then add these fossils to the tree as sampled
  ancestors. These SAs are added as 0-length branches if the `saFormat`
  input is set to `"branch"`, or as degree-2 nodes if it is set to
  `"node"`. The `returnTrueExt` input optionally drops the tip
  representing the true extinction time of a species, and make the last
  sampled fossil of that species the fossil tip, if set to `FALSE`. Note
  that this last functionality required a limited version of the
  `drop.tip` function of the `APE` package to be copied. See
  `?make.phylo` for reference and credits.

## Changes to vignettes

### `overview`

- The Per Capita method estimation was changed to have a more accurate
  estimate.
- A complex example was added to the end to respond to comments from a
  reviewer in the manuscript.
- Examples of new features for version 1.1 (see above) were added.

## Other additions

- `bin.occurrences` allows for post-hoc binning of fossil occurrences,
  so one can take a fossil record including only true times, i.e. an
  output of `sample.clade(..., returnTrue = TRUE)`, and bin it to
  produce the uncertainty in fossil ages ubiquitous in the true fossil
  record.

## Simple fixes

- `sample.clade.R`: added a small bit on help page to explain the
  complication of using `adFun` with extant species.
- `make.phylo.R`: corrected bug in node labels.

# paleobuddy 1.0.0

This is the first release of `paleobuddy`, an R package dedicated to
flexible simulations of diversification, fossil records, and
phylogenetic trees. Below we list current features, and above sections
will be filled as new features and fixes are implemented.

## Main functions

- `rexp.var` generalizes `rexp` to take any function of time as a rate.
  Also allows for a `shape` parameter, in which case it similarly
  generalizes `rweibull`.
- `bd.sim` simulates species diversification with high flexibility in
  allowed speciation and extinction rate scenarios. Produces a `sim`
  object.
- `sample.clade` simulates fossil sampling. Similar flexibility to
  `bd.sim`, though even more so in the case of age-dependent rates.
- `make.phylo` creates a bifurcating phylogenetic tree as a`phylo`
  object (see [APE](https://CRAN.R-project.org/package=ape)) from a
  `sim` object. Can take fossils to be added as length `0` branches.
- `draw.sim` plots a `sim` by drawing species durations and
  relationships, and optionally adding fossils as time points or ranges.

## Secondary functions

- `find.lineages` creates subsets of a `sim` defined as the clades
  descended from one or more species present in the simulation. Needed
  to e.g. generate phylogenetic trees from simulations with more than
  one starting species.
- `make.rate` creates a purely time-dependent (or constant) rate based
  on optional inputs. Used internally to allow for users to define rate
  scenarios easily in `bd.sim`.
- `phylo.to.sim` creates a `sim` object from a `phylo` object, provided
  the user makes choices to solve ambiguities on bifurcating
  phylogenetic trees.
- `var.rate.div` calculates expected diversity for a given
  diversification rate and set of times. Useful for testing `bd.sim` and
  planning rate scenarios.

## S3 classes

- `sim` a class returned by `bd.sim` and used as an input for many
  functions in the package. Formally, it is a named list of vectors
  recording speciation time, extinction time, status (extant or
  extinct), and parent information for each lineage in the simulation.
  It contains the following methods: \*\* `print` gives some quick
  details about number of extant and total species, and the first few
  members of each vector. \*\* `head` and `tail` return the `sim` object
  containing only a given number of species from the beginning and end
  of its vectors, respectively. \*\* `summary` gives quantitative
  details, e.g. quantiles of durations and speciation waiting times.
  \*\* `plot` plots lineages through time (LTT) plots for births,
  deaths, and diversity. \*\* `sim.counts` counts numbers of births,
  deaths, and diversity for some given time `t`. \*\* `is.sim` checks
  the object is a valid `sim` object. Used internally for error
  checking.

## Data

- `temp` temperature data during the Cenozoic. Modified from
  [RPANDA](https://CRAN.R-project.org/package=RPANDA).
- `co2` CO2 data during the Jurassic. Modified from
  [RPANDA](https://CRAN.R-project.org/package=RPANDA).

## Vignettes

- `overview` gives a reasonably in depth look at the main features of
  the package, including examples of workflows using most available rate
  scenarios, and examples of applications.

## Notes

The question of how to structure time came up a lot during development
of the package. Most of the literature in macroevolution and
paleontology considers absolute geological time, i.e. `t = 0` at the
present and `t = 5` five million years ago. It becomes challenging,
however, to visualize and program complex rates going backwards. As
such, the code is structured such that all functions are considered to
go from `0` to the maximum simulation time `tMax`—i.e. the inverse of
absolute geological time. There is only one exception to this rule, the
`adFun` parameter describing age-dependent distribution of fossil
occurrences in `sample.clade`. In any case, all returned objects in the
package are set to follow absolute geological time, so as to conform to
the literature.

## Known issues

- The `integrate` function occasionally fails when given functions that
  vary suddenly and rashly, usually happening in the case of
  environmentally-dependent rates. I have tracked this error down to
  numerical problems in `integrate`, and testing seems to indicate the
  error does not prevent `integrate` from finding the correct result. As
  such this is not currently something I intend to fix, though if issues
  are found that indicate this could be a `paleobuddy` problem, not an
  `integrate` problem, that could change.

- `paleobuddy` is the first package to implement time-dependent
  parameters for Weibull-distributed waiting times. Since the authors
  currently are not aware of an analytical solution to important
  quantities in the BD process in this case, it is challenging to test
  exactly. Simulation tests indicate pretty strongly that the algorithm
  works, however, with one exception—in cases where shape is
  time-dependent and varies dramatically, especially when close to `0`,
  `rexp.var` seems to have a hard time finding the correct waiting time
  distribution. When maintained within levels generally accepted as
  sensible throughout the literature—around 0.8 to 3, say—, and even a
  reasonable amount outside of that range, tests indicate the algorithm
  functions as it should. A `testthat` routine will be implemented in
  the future to formalize these claims, and this issue is one I plan to
  work on soon, especially if users report it as more prevalent than I
  thought.
