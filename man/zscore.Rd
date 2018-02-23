\name{zscore}
\alias{zscore}

\title{ Z-score of relatedness}

\description{
  Compute the Z-score for given relatedness}

\usage{ zscore(app_pairs, m, sd) }

\arguments{
  \item{app_pairs}{ A vector of relatedness values }
  \item{m}{ Reference mean used in Z-score}
  \item{sd}{ Reference standard deviation used in Z-score}
}

\details{
  The Z-score is defined as
  \eqn{Z=\frac{m_{pairs}- m}{sd}}
  where \eqn{m_{pairs}} is the mean of relatedness values.
}

\value{
  Value of Z-score.
}

\author{ Claire Dandine-Roulland }

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Rousset Relationship Matrix
K <- Rousset(x)$Rousset

# Compute Z-score relatedness between individuals from similar population compared to all pairs
zscore(app.pairs(K, x@ped$population), m=app.unpairs(K, x@ped$population, mean), sd=app.unpairs(K, x@ped$population, sd))
}

\seealso{ \code{\link{app.pairs}}, \code{\link{app.unpairs}} }

\keyword{ z-score, Relatedness, pairs }
