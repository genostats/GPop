\name{app.permut}
\alias{app.permut}

\title{ Pairs permuations }

\description{
  Apply given function to relatedness of permutated pairs}

\usage{ app.permut(app, permut.pairs, group=rep(1:nrow(app)), FUN, unpairs=FALSE, fixed.pairs=NULL, B=1e5, thread=1) }

\arguments{
  \item{app}{ Relatedness matrix }
  \item{permut.pairs}{ Vector defining permutated pairs }
  \item{group}{ Vector defining potential group in which permutations are done }
  \item{FUN}{ Function to apply to relatedness of individulas excepted for given pairs. Argument function must be a vector. }
  \item{unpairs}{ Logical, if \code{TRUE} relatedness of individuals pairs excepted for given pairs are considered. }
  \item{fixed.pairs}{ Vector or matrix defining not permutated pairs or a \code{list} of such elements}
  \item{B}{ Number of permutations }
  \item{thread}{ Number of thread used by \code{parallel} package for permutations }
}

\details{
  In each group, pairs definition are permutated in main to apply \code{FUN} function.
  
  It is possible to include fixed rules for paris definition in \code{fixed.pairs} argument.
}

\value{
  A list of two.
  
  First element, \code{app_pairs}, is the value of function \code{FUN} results for relatedness between given pairs or of individulas
  pairs excepted for given pairs for given pairs.
  
  Second element, \code{permut}, is a vector of function \code{FUN} results for relatedness between given pairs or of individulas
  pairs excepted for given pairs for each permutations.
}

\author{ Claire Dandine-Roulland }

\references{Laurent, R and Chaix, R, 2012, \emph{MHC-dependent mate choice in humans: Why genomic patterns from the HapMap European American dataset support the hypothesi}, BioEssays \bold{34(4), 267-271}}

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Rousset Relationship Matrix
K <- Rousset(x)$Rousset

# Permutations
str(app.permut(K, x@ped$population, FUN=mean, B=20, thread=2))
}

\seealso{ \code{\link{app.pairs}}, \code{\link{app.unpairs}}, \code{\link{zscore}} }

\keyword{ Relatedness, pairs, permutations }
