\name{Rousset}
\alias{Rousset}

\title{ Rousset relatedness coefficient Matrix}

\description{
  Compute the Rousset relatedness coefficient Matrix}

\usage{ Rousset(x, which.snps, het_het=0.5, share_het_het=1, autosome.only = TRUE, chunk = 1L) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{which.snps}{ Logical vector, giving which snps to use in the computation.  The default is to use all autosomal SNPs}
  \item{het_het}{ Proportion value of identical variants between pairs of heterozygous individuals}
  \item{share_het_het}{ Number of alleles sharing between pairs of heterozygous individuals}
  \item{autosome.only}{ If \code{TRUE}, only autosomal SNPs will be considered. }
  \item{chunk}{ Parameter for the parallelization: how many SNPs are treated by each task}
}

\details{
  The Rousset relatedness coefficient between two individuals i and j, is defined as
  \eqn{R_{ij}=\frac{Q_{ij}- Q_m}{1-Q_m}}
  where \eqn{Q_{ij}} is the proportion of variants which are identical by state between the two individuals and \eqn{Q_m} is the
  mean proportion of variants which are identical by state in the population. At a given SNP, the proportion of identical variants
  for a given pair of individuals was set to 0 if both individuals were homozygous and carrying a different allele (e.g. 00 and 11),
  1 if both individuals were homozygous and carrying the same allele (e.g. 00 and 00), het_het for pairs of heterozygous individuals
  (e.g. 01 and 01) and 0.5 in all others cases.
}

\value{
  A list of seven.
  
  \code{Rousset} element is the Rousset relatedness coefficient Matrix which is a symmetric square matrix of dimension equal to the number of individuals.
  
  \code{F} element contains the matrix of \eqn{Q_{ij}} values given in details section.
  
  \code{Qm} element is the \eqn{Q_m} value given in details section.
  
  \code{Nsnp} element is a symmetric square matrix of dimension equal to the number of individuals which contains the number of SNPs available
              by pairs
  
  \code{SHare0} element is the Matrix of SNP proportion sharing 0 allele between two individuals which is a symmetric square matrix of dimension
                equal to the number of individuals
  
  \code{SHare1} element is the Matrix of SNP proportion sharing 1 allele between two individuals which is a symmetric square matrix of dimension
                equal to the number of individuals. By default, pairs of heterozygous individuals are considered to share 1 allele.
  
  \code{SHare2} element is the Matrix of SNP proportion sharing 2 alleles between two individuals which is a symmetric square matrix of dimension
                equal to the number of individuals
}

\author{ Claire Dandine-Roulland }

\references{Rousset, F, 2002, \emph{Inbreeding and relatedness coefficients: what do they measure?}, Heredity \bold{88(5), 371-380}}

\examples{
# load chr2 data set (~10k SNPs in low LD)
x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )

# Compute Rousset Relationship Matrix
K <- Rousset(x)
str(K)
dim(K)
}

\keyword{ Rousset relatedness coefficient }
