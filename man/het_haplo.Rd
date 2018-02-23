\name{het_haplo}
\alias{het_haplo}

\title{ Haplotype heterozygoty }

\description{ Global and by SNP Fst computation for given populations }
\usage{ het_haplo(x, map=c('B36','B37'), block=0.5, minl=1, maxl=+Inf, thread)}

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
het <- het_haplo(x, map='B37', block=0.5, minl=1, maxl=+Inf, thread
str(het)
}

\keyword{ FST }
