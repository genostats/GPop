\name{Fst}
\alias{Fst}

\title{ Fst Computation }

\description{ Global and by SNP Fst computation for given populations }
\usage{ Fst( x, pop, beg = 1, end = ncol(x) )}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{pop}{ The population vector. }
  \item{beg}{ Index of the first SNP used for Fst computation }
  \item{end}{ Index of the last SNP used for Fst computation }
}

\value{
  A List of two.
  
  First element is a data frame, giving for each considered SNP, some of the following columns :
  \item{a, b and c} {Fst computation elements defined by Weir and Cockerham}
  \item{Fst} {Fst}
  Second element gives global Fst computed on all considered SNPs.
}

\author{ Claire Dandine-Roulland }

\references{Weir BS and Cockerham CC, 1984, \emph{Estimating F-statistics for the analysis of population structure}, Evolution \bold{38(6), 1358-1370}}

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- set.stats(read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") ))

# Compute Fst
fst <- Fst(x, x@ped$population)
str(fst)
}

\keyword{ FST }