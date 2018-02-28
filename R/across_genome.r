app_genome_test <- function(bed, FUN=mean, relatedness=c('GRM', 'Rousset'), length=3600000, wind.bp=300000, length.snp=NULL, wind.snp=NULL, thread=22, mapgen='All_Hapmap', build='B36', LD=NULL, het_het=0.5)
{
  # load of genetic map
  if (build=='B36') {
    if (mapgen %in% c('All_Hapmap', 'CEU_YRI_Hapmap', 'CEU_1000G') ) map <- genmapB36[[mapgen]] else stop("'mapgen' option is not valid")
  } else if (build=='B37') {
    if (mapgen %in% c('All_1000G') ) map <- genmapB37[[mapgen]] else stop("'mapgen' option is not valid")
  } else stop("'build' option is not valid")

  # load centromere place for each chromosom
  centro <- centromere[[build]]
  
  if (thread>1) { RNGkind("L'Ecuyer-CMRG"); cl <- makeForkCluster(nnodes=thread); }
  
  if (is.character(bed))
  {
    
    snps <- read.table(paste(bed, '.bim', sep=''))
  } else if (is.bed.matrix(bed))
  {
    snps <- bed@snps
  } else {
    stopCluster(cl)
	stop('"bed" must be a string or a bed.matrix.') 
  }
 
  result <- NULL
  for (i in snps$chr)
  {
    if (is.null(wind.bp))
	{
	  w <- snps[snps$chr==i,]
	  temp <- data.frame(chr=i, start=w$pos[seq(1, floor(nrow(w)/wind.snp)*wind.snp, by=wind.snp)] )
      temp$end <- w$pos[seq(1, floor(nrow(w)/wind.snp)*wind.snp, by=wind.snp)+length.snp]
      temp$end[length(temp$end)] <- w$pos[nrow(w)]
	  temp$index.start <- seq(1, floor(nrow(w)/wind.snp)*wind.snp, by=wind.snp)
	  temp$index.end <- seq(1, floor(nrow(w)/wind.snp)*wind.snp, by=wind.snp)+length.snp
	  temp$recombi_all <- apply(cbind(temp$start, temp$end), 1, function(x) mean( map[[paste('chr', i, sep='')]]$rate.cM.Mb[ map[[paste('chr', i, sep='')]]$base>=x[1] & map[[paste('chr', i, sep='')]]$base<=x[2] ] ) )
	  temp$recombi_region <- apply(cbind(temp$start, temp$end), 1, function(x) {
        n <- which(map[[paste('chr', i, sep='')]]$base<=x[2] & map[[paste('chr', i, sep='')]]$base>=x[1]); ifelse(
	                         length(n)==1, map[[paste('chr', i, sep='')]]$rate.cM.Mb[n],
							 ( map[[paste('chr', i, sep='')]]$cM[ max(w) ] - map[[paste('chr', i, sep='')]]$cM[ w[1] ] )/( map[[paste('chr', i, sep='')]]$base[ max(w) ]/1e6 - map[[paste('chr', i, sep='')]]$base[ w[1] ]/1e6 ));})
    } else {
	  w <- snps[snps$chr==i,]
      temp <- data.frame(chr=i, start=seq(0, as.integer(ceiling( max(w$pos)/wind.bp )*wind.bp), by=wind.bp) )
      temp$end <- temp$start + length
 	  temp$index.start <- sapply(seq(0, as.integer(ceiling( max(w$pos)/wind.bp )*wind.bp), by=wind.bp), function(x) which.min(ifelse(w$pos>=result$start[i], w$pos, Inf)) )
	  temp$index.end <- sapply(seq(0, as.integer(ceiling( max(w$pos)/wind.bp )*wind.bp), by=wind.bp), function(x) which.max(ifelse(w$pos<result$end[i], w$pos, -Inf)) )
	  temp$recombi_all <- apply(cbind(temp$start, temp$end), 1, function(x) mean( map[[paste('chr', i, sep='')]]$rate.cM.Mb[ map[[paste('chr', i, sep='')]]$base>=x[1] & map[[paste('chr', i, sep='')]]$base<=x[2] ] ) )
	  temp$recombi_region <- apply(temp[,2:3], 1, function(x) {
        n <- which(map[[paste('chr', i, sep='')]]$base<=x[2] & map[[paste('chr', i, sep='')]]$base>=x[1]); ifelse(
	                         length(n)==1, map[[paste('chr', i, sep='')]]$rate.cM.Mb[n],
							 ( map[[paste('chr', i, sep='')]]$cM[ max(n) ] - map[[paste('chr', i, sep='')]]$cM[ n[1] ] )/( map[[paste('chr', i, sep='')]]$base[ max(n) ]/1e6 - map[[paste('chr', i, sep='')]]$base[ n[1] ]/1e6 ));})
   }
    result <- rbind(result, temp)
  }
  rm(temp,x,w)
  result$centro <- ( (result$start > centro$start[result$chr] & result$start < centro$end[result$chr]) | (result$end > centro$start[result$chr] & result$end < centro$end[result$chr]) | (result$end > centro$end[result$chr] & result$start < centro$start[result$chr]) )
 
  #if(!file.exists(paste(bed, 'bed',sep='.')))
  
  if (thread>1)
  {
    e <- environment()
	
	withGlobals <- function(FUN, ...){ 
    environment(FUN) <- list2env(list(...)) 
    FUN 
} 

    # create cluster
    RNGkind("L'Ecuyer-CMRG")
    cl <- makeCluster(thread)
    clusterSetRNGStream(cl, runif(1)*(2**31-1) )
    clusterCall(cl, function() library(gaston))
    clusterCall(cl, function() library(gaston.pop))
    clusterExport(cl, c("map", "LD", "FUN", "e", "result"), envir=environment())

	r <- parRapply(cl, cbind(result$chr, result$start, result$end), withGlobals( function(i) {
	  #x <- read.bed.matrix.part(bed, beg=result$index.start[i], end=result$index.end[i])
	  print('ici')
	  t <- list()
      t$num <- ncol(x)
      if (ncol(x)>0) {
        standardize(x) <- 'p'
        if (!is.null(LD) & ncol(x)>1) l <- LD.thin(x, LD, extract=TRUE)
        t$num <- ncol(x)
        #if (sum(map[[paste('chr', unique(x@snps$chr), sep='')]]$base %in% x@snps$pos)>0)  t$recombi <- map[[paste('chr', unique(x@snps$chr), sep='')]]( map[[paste('chr', unique(x@snps$chr), sep='')]]$rate.cM.Mb[ map[[paste('chr', unique(x@snps$chr), sep='')]]$base %in% x@snps$pos ], na.rm=TRUE)
      }

      # LD
      L <- LD(x, lim=c(1,ncol(l)))
	  L <- L[lower.tri(L, diag=FALSE)]
      t$LD <- mean(L)
      t$LD_sd <- sd(L)
      t$maf <- mean(x@snps$maf)
      t$maf_sd <- sd(x@snps$maf)
     
      if ('GRM' %in% relatedness)
      {
        app <- GRM(x)
		for (f in FUN) t[[paste('grm', deparse(quote(f)), sep='_')]] <- f(app, x)
      }
      if ('Rousset' %in% relatedness)
      {
        app <- Rousset(x, het_het)$Rousset
		for (f in FUN) t[[paste('rousset', deparse(quote(f)), sep='_')]] <- f(app, x)
      }
      return(result) }, x=bed) )
    stopCluster(cl)
  } else r <- NA
  
  return(r)
}
