\name{app.unpairs}
\alias{app.unpairs}

\title{ Relatedness excepted given pairs }

\description{
  Removed given pairs from relatedness matrix and apply given function}

\usage{ app.unpairs(app, pairs, FUN=identity) }

\arguments{
  \item{app}{ Relatedness matrix }
  \item{unpairs}{ Vector or matrix defining wanted pairs or a \code{list} of such elements}
  \item{FUN}{ Function to apply to relatedness of individulas excepted for given pairs. Argument function must be a vector. }
}

\value{
  Result value of FUN function for relatedness of individuals pairs excepted for given pairs.
}

\author{ Claire Dandine-Roulland }

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Rousset Relationship Matrix
K <- Rousset(x)$Rousset

# Cmmpute mean of relatedness between individuals from different population
app.unpairs(K, x@ped$population, FUN=mean)
}

\seealso{ \code{\link{app.pairs}} }

\keyword{ Relatedness, pairs }
