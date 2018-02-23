\name{het.haplo}
\alias{het.haplo}

\title{ Haplotype heterozygosity rate }

\description{ Global and by SNP Fst computation for given populations }
\usage{ het.haplo(x, map=c('B36','B37'), block=0.5, minl=1, maxl=+Inf, thread=1)}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{map}{ Recombinaison map. Must be \code{B36}, \code{B37} or a data.frame contains 'chr', 'base', 'rate.cM.Mb' and 'cM' variables".
              If \code{B36}, we use recombinaison map from Hapmap contained in \code{genmapB36}.
			  If \code{B37}, we use recombinaison map from 1000 Genomes contained in \code{genmapB37}.}
  \item{block}{ Maximum recombinaison value in cM/Mb to define haplotype blocks }
  \item{minl}{ Minimum length of haplotype block }
  \item{maxl}{ Maximum length of haplotype block }
  \item{thread}{ Number of used thread }
}

\value{
  A List of 22 elements, one by autosomal chromosoms.
  Each element is a list of two with:
  - a vector with heterozygotie rate by considered haplotype blocks.
  - a list which gives for each haplotype blocks ids of considered SNP.
}

\details{
  Haplotype frequencies are estimated using EM algorithm available in \code{\link{haplo.stats}} package.
}

\author{ Claire Dandine-Roulland }

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- set.stats(read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") ))

# Compute Fst
het <- het.haplo(x, map='B37', block=0.5, minl=1, maxl=+Inf, thread=1)
str(het)
}

\keyword{ Haplotype Heterozygosity rate }
