Fst <- function(x, pop, beg = 1, end = ncol(x))
{  
  if(beg < 1 || end > ncol(x)) stop("range too wide")
  if(is.null(x@snps$hz)) stop("Need hz to be set in x")
  if(is.null(x@p)) stop("Need p to be set in x")
  if (nrow(x)!=length(pop)) stop("Dimensions of pop and x mismatch")
  
  # heterozygotie and frequencies on all populations
  pb <- x@p[beg:end]
  hb <- x@snps$hz[beg:end]
  
  # Number of populations
  r <- length(unique(pop))
  
  # c (for each SNP)
  c <- hb/2
  
  # Compute genotype counts in each population
  N <- .Call("gg_Fst_countbypop", PACKAGE = "gaston.pop", x@bed, as.numeric(as.factor(pop)), beg-1, end-1)
  
  # n
  ni <- N$N0+N$N1+N$N2
  nc <- ( rowSums(ni) - rowSums(ni**2)/rowSums(ni) )/(r-1)
  nb <- rowSums(ni)/r
  
  # heterozygotie and frequencies in each populations
  p <- (N$N1+2*N$N2)/(2*ni)
  h <- N$N1/ni
  
  # Compute weighted variance of allele frequencies
  s2 <- rowSums( ni*(p-pb)**2 )/( (r-1)*nb )
    
  # a (for each SNP)
  a <-  nb/nc * ( s2 - ( pb*(1-pb) - (r-1)/r*s2 - hb/4 )/(nb-1) )

  # b (for each SNP)
  b <-  nb/(nb-1) * ( pb*(1-pb) - (r-1)/r*s2 - (2*nb-1)/(4*nb)*hb )
  
  # SNP information
  L <- list(chr = x@snps$chr, pos = x@snps$pos, id  = x@snps$id)
  if(beg > 1 | end < ncol(x))  # avoid copy
  L <- L[beg:end,] 
 
  return(list(bySNP=data.frame( L, a=a, b=b, c=c, fst=a/(a+b+c)), global=sum(a)/sum(a+b+c)))
}