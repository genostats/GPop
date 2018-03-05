\name{app.genome}
\alias{app.genome}

\title{ Across the genome }

\description{
  Apply given function to relatedness at sliding genomic windows }

\usage{ app.genome(bed, FUN, relatedness=GRM, windows, sliding,
                   unit=c("bases", "indices"), map='B36', centromere=NULL, LD.thin=NULL, thread=22) }

\arguments{
  \item{bed}{ A \code{bed.matrix}}
  \item{FUN}{ Function to apply to relatedness of individulas at genomic windows. Argument function must be a relatedness matrix. }
  \item{relatedness}{ Function to compute relatedness matrix at genomic windows. Argument function must be a bed.matrix. }
  \item{windows}{ Numerical value defined length of genomic windows }
  \item{sliding}{ Increment length value of sliding genomic windows }
  \item{unit}{ Distance unit in \code{windows} and \code{sliding} }
  \item{map}{  Recombinaison map. Must be \code{B36}, \code{B37} or a data.frame contains 'chr', 'base', 'rate.cM.Mb' and 'cM' variables".
              If \code{B36}, we use recombinaison map from Hapmap contained in \code{genmapB36}.
			  If \code{B37}, we use recombinaison map from 1000 Genomes contained in \code{genmapB37}.
			  If \code{data.frame}, \code{centromere} argument must be given.}
  \item{centromere}{ Optional \code{data.frame} contains 'chr', 'start', 'end' variables which define centromere localisations (one line by chromosom) }
  \item{LD.thin}{ Distance unit in \code{windows} and \code{sliding} }
  \item{thread}{ Number of thread used by \code{parallel} package for permutations }
}

\details{
  Compute the genetic relatedness and apply given function(s) for all sliding windows of given length \code{windows} in increments of \code{sliding}.

\value{
  A data.frame (one line by genomic region) with columns :
  
  \code{chr}, \code{start} end \code{end} defining genomic region localisation

  \code{index.start} end \code{index.end} defining first and last SNP index inluded in genomic region

  \code{recombi_snps}, recombinaison rate estimated on SNPs included in \code{x} defining genomic region

  \code{recombi_region}, mean recombinaison rate on genomic region

  \code{centro}, \code{logical}, \code{TRUE} if region overlapping a centromere

  \code{num}, number of SNPs in genomic region

  \code{LD} end \code{LD_sd}, mean and standard deviation of LD in genomic region
  
  \code{maf} end \code{maf_sd}, mean and standard deviation of maf in genomic region
}

\author{ Claire Dandine-Roulland }

\references{Laurent, R and Chaix, R, 2012, \emph{MHC-dependent mate choice in humans: Why genomic patterns from the HapMap European American dataset support the hypothesi}, BioEssays \bold{34(4), 267-271}}

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Rousset Relationship Matrix
relatedness <- function(x) Rousset(x)$Rousset
fun <- function(app) app.pairs(app, x@ped$population, mean)

# Across the genome
head(app.genome(x, fun, relatedness, windows=2000, sliding=500, unit="bases", thread=2))
}

\seealso{ \code{\link{app.pairs}}, \code{\link{app.unpairs}}, \code{\link{zscore}} }

\keyword{ Relatedness, genome, sliding windows }
