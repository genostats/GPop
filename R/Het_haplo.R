# Extract variables needed in haplo.stats
bedtogeno <- function(bed)
{
  standardize(bed) <- 'none'
  locus <- bed@snps$id
  bed <- as.matrix(bed)
  geno <- matrix(NA, nrow=nrow(bed), ncol=2*ncol(bed))

  for (i in 1:nrow(bed)) for (j in 1:ncol(bed)) if (!is.na(bed[i,j])) geno[i, ((j-1)*2+1):(2*j)] <- c(1*(bed[i,j] %in% c(1,2)), 1*(bed[i,j]==2))
  
  return(list(geno=geno, locus=locus))
}


#### Haplotype Heterozygosity 

het_haplo <- function(x, map=c('B36','B37'), block=0.5, minl=1, maxl=+Inf, thread)
{
  # load of genetic map
  if (mapgen %in% c('B36','B37') ) { data(paste('genmap', map, sep='')); map <- eval(paste(text=paste('genmap', map, sep='')))$genmap; }
  else if (is.data.frame(map)) if (names(map)!=c("chr", "base", "rate.cM.Mb", "cM")) stop("'map' must contains 'chr', 'base', 'rate.cM.Mb' and 'cM' variables")
  else stop("'map' must be equal to 'B36' or 'B37' or a data frame")

  if (thread>1)
  {
    if (thread==2) xx <- c(1,22,2,21,3,20,4,19,5,18,6,17,7,16,8,15,9,14,10,13,11,12)
    else if (thread==3) xx <- c(1,22,21,20,2,19,18,3,17,16,4,15,14,5,13,12,6,11,10,7,9,8)
    else if (thread==5) xx <- c(1,22,21,20,19,2,18,17,16,15,3,14,13,12,4,11,10,9,5,8,7,6)
    else if (thread==7) xx <- c(1,22,21,20,2,19,18,3,17,16,4,15,14,5,13,12,6,11,10,7,9,8)
    else if (thread==11) xx <- c(1,22,2,21,3,20,4,19,5,18,6,17,7,16,8,15,9,14,10,13,11,12)
    else xx <- 1:22

    # create cluster
    RNGkind("L'Ecuyer-CMRG")
    #if (!is.character(bed))
    #{
      cl <- makeForkCluster(nnodes=thread)
      clusterSetRNGStream(cl, runif(1)*(2**31-1) )
    #} else {
    #  cl <- makeCluster(thread)
    #  clusterSetRNGStream(cl, runif(1)*(2**31-1) )
    #  clusterCall(cl, function() library(gaston))
    #  clusterCall(cl, function() library(gaston.pop))
    #  clusterCall(cl, function() { library(RhpcBLASctl); blas_set_num_threads(1); })
    #  clusterExport(cl, c("app_genome_chr", "app_couple", "z_score", "mean_couple"), envir=globalenv())
    #  clusterExport(cl, c("ind", "snps", "map", "relatedness", "method", "bed", "length", "wind.bp", "length.snp", "wind.snp", "centro", "LD"), envir=environment())
    #}
    
    r <- parLapply(cl, xx, function(i) {
      if (is.character(bed))
      {
        if (file.exist(paste(bed, i, sep=''))) chr <- read.bed.matrix(paste(bed, i, sep='')) else {
          warnings(paste('File ', bed, i,  "doesn't exist.", sep=''))
          return(NULL)}
      } else chr <- select.snps(bed, chr==i)
      if (ncol(chr)==0) {warnings(paste('None SNP i chromosome', i)); return(NULL);}
      t <- map[map$chr==i,2:4]
      het_haplo_chr(chr, t, block=block, minl=minl, maxl=maxl, method=method, plink=plink)} )
    stopCluster(cl)
    r <- r[order(xx)]
  } else {
    r <- list()
    for (i in 1:22)
    {
      if (is.character(bed))
      {
        if (file.exist(paste(bed, i, sep=''))) chr <- read.bed.matrix(paste(bed, i, sep='')) else {
          warnings(paste('File ', bed, i,  "doesn't exist.", sep=''))
          next }
      } else chr <- select.snps(bed, chr==i)
      if (ncol(chr)==0) {warnings(paste('None SNP i chromosome', i)); next;}
      print(i)
      t <- map[map$chr==i,2:4]
      r[[i]] <- het_haplo_chr(chr, t, block=block, minl=minl, maxl=maxl, method=method, plink=plink)
    }
  }
}
         

## By chromosome
het_haplo_chr <- function(chr, map, block=0.5, minl=1, maxl=+Inf)
{
  # Verify order of SNP
  if (is.unsorted(chr@snps$pos)) chr <- chr[, order(chr@snps$pos)]

  # Selection of SNPs commun in chr and map
  temp <- chr@snps$pos %in% map$base
  tab <- data.frame(snps=chr@snps$id[temp], pos=chr@snps$pos[temp], cM_Mb=NA)

  # take mean recombinaison rate betwwen known SNPs
  temp <- match(tab$pos, matrix(map$base, ncol=1))
  temp <- temp[temp!=F]
  for (i in 1:(nrow(tab)-1)) tab$cM_Mb[i+1] <- mean( map$rate.cM.Mb[ (temp[i]+1):temp[i+1] ] )

  # Define blocks
  temp <- which( tab$cM_Mb<block ) # limit to define block
  j <- 1
  k <- 1
  b <- list()
  while (j<length(temp))
  {
    #b[[k]] <- as.character(tab$snps[j])
    b[[k]] <- NA
    i <- j
    while( temp[i+1]==(temp[i]+1) & i<length(temp) ) {b[[k]] <- c(b[[k]], as.character(tab$snps[temp[i]])); i<-i+1;}
    k <- k+1
    j <- i+1
  }
  
  # create plink file for block
  b <- lapply(b, function(x) x[-1])
  b <- b[which(unlist(lapply(b, length))>=minl & unlist(lapply(b, length))<=maxl)]

  het <- rep(0, length(b))
  for (i in 1:length(b))
  {
    temp <- select.snps(chr, id %in% b[[i]])
    temp <- bedtogeno(temp)
    het[i] <- 2*nrow(chr)/(2*nrow(chr)-1)*(1-sum(haplo.em(temp$geno, temp$locus, miss.val=NA)$hap.prob**2))
  }
  return(list(het=het, blocks=b))
}




