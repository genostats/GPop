Rousset <- function(x, which.snps, het_het=0.5, autosome.only = TRUE, chunk = 1L) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  if(!is.numeric(het_het) | het_het<0 | het_het>1)
    stop("het_het must be between 0 and 1")

  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")

  K <- .Call('gg_Rousset', PACKAGE = 'gaston.pop', x@bed, which.snps, het_het, chunk) 

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0) 
      rownames(K$Rousset) <- colnames(K$Rousset) <- x@ped$id
      rownames(K$F) <- colnames(K$F) <- x@ped$id
  }

  K
}
