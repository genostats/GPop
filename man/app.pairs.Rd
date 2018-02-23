\name{app.pairs}
\alias{app.pairs}

\title{ Relatedness between pairs }

\description{
  Extract given pairs from relatedness matrix and apply given function}

\usage{ app.pairs(app, pairs, FUN=identity) }

\arguments{
  \item{app}{ Relatedness matrix }
  \item{pairs}{ Vector or matrix defining wanted pairs or a \code{list} of such elements}
  \item{FUN}{ Function to apply to relatedness of given pairs. Argument function must be a vector. }
}

\value{
  Result value of FUN function for relatedness between given pairs.
}

\author{ Claire Dandine-Roulland }

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Rousset Relationship Matrix
K <- Rousset(x)$Rousset

# Compute mean of relatedness between individuals from similar population
app.pairs(K, x@ped$population, FUN=mean)
}

\seealso{ \code{\link{app.unpairs}} }

\keyword{ Relatedness, pairs }
