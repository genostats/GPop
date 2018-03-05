nb.snps <- function(basename) {
  bim <- paste(basename, ".bim", sep="")
  rds <- paste(basename, ".rds", sep="")
  if(!file.exists(bim) & !file.exists(rds)) { # peut-être on a donné le .bed pour basename
    if(grep("\\.bed$", basename))
      basename <- sub("\\.bed$", "", basename)
    bim <- paste(basename, ".bim", sep="")
    rds <- paste(basename, ".rds", sep="")
  }
  if(!file.exists(bim) & !file.exists(rds)) 
    stop("No file ", rds, " or ", bim)

  if(!file.exists(rds)) 
    nb.snps <- R.utils::countLines(bim)
  else {
    x <- readRDS(rds)
    nb.snps <- nrow(x@snps) 
  }
  nb.snps
}

read.bed.matrix.part <- function(basename, bed = paste(basename, ".bed", sep=""), fam = paste(basename, ".fam", sep=""), 
                                 bim = paste(basename, ".bim", sep=""), rds = paste(basename, ".rds", sep=""), beg = 1, end,
                                 verbose = getOption("gaston.verbose",TRUE)) {

  bed <- path.expand(bed)
  if(!file.exists(bed)) { # peut-être on a donné le .bed pour basename
    if(grep("\\.bed$", basename)) {
      basename <- sub("\\.bed$", "", basename)
      bim <- paste(basename, ".bim", sep="")
      fam <- paste(basename, ".fam", sep="")
      bed <- paste(basename, ".bed", sep="") 
      if(!is.null(rds)) rds <- paste(basename, ".rds", sep="")
    }
  }

  if(!file.exists(bed)) stop("file ", bed, " not found")

  if(is.null(rds) || !file.exists(rds)) {
    if(!file.exists(fam)) stop("file ", fam, " not found")
    if(!file.exists(bim)) stop("file ", bim, " not found")

    if(verbose) cat("Reading", fam, "\n")
    ped <- read.table(fam, stringsAsFactors = FALSE) 
    colnames(ped) <- gaston:::pednames

    if(verbose) cat("Reading", bim, "\n")
    snp <- read.table(bim, stringsAsFactors = FALSE)
    colnames(snp) <- gaston:::snpnames

    if(missing(end)) 
      end <- nrow(bim)

    if(beg > end | beg > nrow(snp) | beg < 1 | end > nrow(snp)) 
      stop("beg (", beg, ") and end (", end, ") values incompatible with ", bim, " file (", nrow(bim), " rows)")

    if(verbose) cat("Reading", bed, "\n")
    bed <- .Call('read_bed_file_part', PACKAGE = "gaston.pop", bed, nrow(ped), beg, end)

    snp <- snp[beg:end, ]
    x <- new("bed.matrix", bed = bed, snps = snp, ped = ped,                            
      p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
      standardize_mu_sigma = FALSE )
    if(getOption("gaston.auto.set.stats", TRUE)) x <- gaston::set.stats(x, verbose = verbose)
    return(x)
  }
 
  # reading rds [on ne fait pas set.stats dans ce cas !!]
  if(verbose) cat("Reading", rds, "\n")
  x <- readRDS(rds)
  if ( is(x) != "bed.matrix" ) stop("The object in file ", rds, " is not a bed.matrix")

  if(missing(end)) 
    end <- nrow(x@snps)
  
   if(beg > end | beg > nrow(x@snps) | beg < 1 | end > nrow(x@snps))
      stop("beg (", beg, ") and end (", end, ") values incompatible with ", rds, "file (", nrow(x@snps), " snps)")


  if(verbose) cat("Reading", bed, "\n")
  x@bed <- .Call('read_bed_file_part', PACKAGE = "gaston.utils", bed, nrow(x@ped), beg, end)
  x@snps <- x@snps[beg:end,]
  x@p  <- x@p[beg:end]
  x@mu <- x@mu[beg:end]
  x@sigma <- x@sigma[beg:end]
  x
}


