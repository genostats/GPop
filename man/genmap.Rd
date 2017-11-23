\name{genmap}
\docType{data}
\alias{genmapB36}
\alias{genmapB37}

\title{Genetic map data set}

\description{
  The data set contains genetic maps from Hapmap or 1000 Genomes projects for build 36 and build 37 or GRCh37 respectively.
}
\usage{ data(genmapB36); data(genmapB37); }

\format{
  There are two lists (\code{genmapB36} for NCBI36 assembly and \code{genmapB37} for GRCh37 assembly) of two:
  \describe{
    \item{\code{genmap}}{ Data frames containing chromosom, base position, cM/Mb rate and cM distance by chromosome.}
     \item{\code{centromere}}{ For each chromosome the beguinning and the end of centromere }
  }
}

\source{ The data were obtained from Hapmap project (see \url{ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/}) and the 1000 Genomes project
         (see \url{ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/}). }

\examples{
data(genmap)
str(genmap)
str(centromere)
}

\keyword{Genetic maps}