Rousset <- function(x, which.snps, het_het=0.5, share_het_het=1, autosome.only = TRUE, chunk = 1L) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  if(!is.numeric(het_het) | het_het<0 | het_het>1)
    stop("het_het must be between 0 and 1")
    
  if(!( share_het_het %in% 0:2))
    stop("share_het_het must be egal to 0, 1 or 2")
  if (share_het_het==0) {share0_het_het <- 1; share1_het_het <- share2_het_het <- 0;}
  else if (share_het_het==1) {share1_het_het <- 1; share0_het_het <- share2_het_het <- 0;}
  else if (share_het_het==2) {share2_het_het <- 1; share0_het_het <- share1_het_het <- 0;}

  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")

  K <- .Call('gg_Rousset', PACKAGE = 'GPop', x@bed, which.snps, het_het, share0_het_het, share1_het_het, share2_het_het, chunk) 

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0) 
      rownames(K$Rousset) <- colnames(K$Rousset) <- x@ped$id
      rownames(K$F) <- colnames(K$F) <- x@ped$id
      rownames(K$Nsnp) <- colnames(K$Nsnp) <- x@ped$id
      rownames(K$Share0) <- colnames(K$Share0) <- x@ped$id
      rownames(K$Share1) <- colnames(K$Share1) <- x@ped$id
      rownames(K$Share2) <- colnames(K$Share2) <- x@ped$id
  }

  K
}

#reshape. <- function(K, include = c(-Inf, +Inf), exclude) {
#  diag(K) <- NA
#  if(missing(exclude))
#    w <- which(include[1] < K & K < include[2])
#  else 
#    w <- which(include[1] < K & K < include[2] & (K < exclude[1] | K > exclude[2]))
#  I <- row(K)[w]
#  J <- col(K)[w]
#  R <- K[w]
#  ww <- (I < J)
#  i <- I[ww];
#  j <- J[ww];
#  data.frame(i = i, j = j, id_i = rownames(K)[i], id_j = colnames(K)[j], k = R[ww])
#}

