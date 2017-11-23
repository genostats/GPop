\name{genmap}
\docType{data}
\alias{genmap}
\alias{centromere}

\title{Genetic map data set}

\description{
  The data set contains genetic maps from Hapmap or 1000 Genomes projects for build 36 and build 37 or GRCh37 respectively.
}
\usage{ data(genmap) }

\format{
  There are two data objects in the dataset:
  \describe{
    \item{\code{genmap}}{ List of genetic maps for NCBI36 and GRCh37 assemblies as a data frames containing chromosom, base position, cM/Mb rate
	                         and cM distance by chromosome.}
     \item{\code{centromere}}{ List of two data frames containing for each chromosome the beguinning and the end of centromere for
                              build 36 (\code{$B36}) and build 37 or GRCh37 (\code{$B37})}
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