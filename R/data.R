#' Cenozoic temperature data
#' 
#' Temperature data during the Cenozoic. Modified from the \code{InfTemp} data
#' set in \href{https://github.com/hmorlon/PANDA}{RPANDA}, originally inferred 
#' from delta O18 measurements.
#' 
#' @usage 
#' 
#' data(temp)
#' 
#' @format A data frame with 17632 rows and 2 variables:
#' \describe{
#'   \item{t}{A numeric vector representing time since the beginning of the data 
#'   frame age, approximately 67 million years ago, in million years. We set this
#'   from past to present as opposed to present to past since birth-death 
#'   functions in \code{paleobuddy} consider time going in the former direction.}
#'   
#'   \item{temperature}{A numeric vector representing temperature in degrees 
#'   celsius corresponding to time \code{t}. Note there might be more than one 
#'   temperature for each time \code{t} given the resolution of the data set.}
#' }
#' 
#' @source \url{https://github.com/hmorlon/PANDA}
#' 
#' @references 
#' 
#' Morlon H. et al (2016) RPANDA: an R package for macroevolutionary analyses on 
#' phylogenetic trees. \emph{Methods in Ecology and Evolution} 7: 589-597.
#' 
#' Epstein, S. et al (1953) Revised carbonate-water isotopic temperature scale 
#' \emph{Geol. Soc. Am. Bull.} 64: 1315-1326.
#' 
#' Zachos, J.C. et al (2008) An early Cenozoic perspective on greenhouse warming 
#' and carbon-cycle dynamics \emph{Nature} 451: 279-283.
#' 
#' Condamine, F.L. et al (2013) Macroevolutionary perspectives to environmental 
#' change \emph{Eco Lett.} 16: 72-85.
"temp"

#' Jurassic CO2 data
#' 
#' CO2 data during the Jurassic. Modified from the \code{co2} set in 
#' \href{https://github.com/hmorlon/PANDA}{RPANDA}, originally taken from Mayhew 
#' et al (2008, 2012). 
#' 
#' @usage 
#' 
#' data(co2)
#' 
#' @format A data frame with 53 rows and 2 variables:
#' \describe{
#'   \item{t}{A numeric vector representing time since the beginning of the data 
#'   frame age, 520 million years ago, in million years. We set this from past to 
#'   present as opposed to present to past since birth-death functions in 
#'   \code{paleobuddy} consider time going in the former direction.}
#'   
#'   \item{co2}{A numeric vector representing CO2 concentration as the ratio of
#'   CO2 mass at \code{t} over the present.}
#' }
#' 
#' @source \url{https://github.com/hmorlon/PANDA}
#' 
#' @references 
#' 
#' Morlon H. et al (2016) RPANDA: an R package for macroevolutionary analyses on 
#' phylogenetic trees. \emph{Methods in Ecology and Evolution} 7: 589-597.
#' 
#' Mayhew, P.J. et al (2008) A long-term association between global temperature
#' and biodiversity, origination and extinction in the fossil record 
#' \emph{Proc. of the Royal Soc. B} 275:47-53.
#' 
#' Mayhew, P.J. et al (2012) Biodiversity tracks temperature over time 
#' \emph{Proc. of the Nat. Ac. of Sci. of the USA} 109:15141-15145.
#' 
#' Berner R.A. & Kothavala, Z. (2001) GEOCARB III: A revised model of atmospheric 
#' CO2 over Phanerozoic time \emph{Am. J. Sci.} 301:182–204.
"co2"
