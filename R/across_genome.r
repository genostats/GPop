app.genome <- function(bed, FUN, relatedness, windows, sliding, unit=c("bases", "indices"), map='B37', centromere=NULL, LD.thin=NULL, thread=22)
{
  ## Verif
  if (!is.function(FUN) & !is.list(FUN)) stop('"FUN" argument must be a function of a list of functions')
  if (is.list(FUN)) for (i in 1:length(FUN)) if (!is.function(FUN[[i]])) stop('"FUN" argument must be a function of a list of functions')
  if (is.function(FUN)) { f <- as.character(substitute(FUN)); FUN <- list(FUN); names(FUN) <- ifelse(length(f)>1, 'function1', f); }
  
  if (!is.function(relatedness)) stop('"relatedness" argument must be a function')
  if (is.list(relatedness)) for (i in 1:length(relatedness)) if (!is.function(relatedness[[i]])) stop(relatedness)
  if (is.function(relatedness)) { f <- as.character(substitute(relatedness)); relatedness <- list(relatedness); names(relatedness) <- ifelse(length(f)>1, 'relatedness1', f); }

  # load of genetic map
  if (map=='B36') {
    data('genmapB36')
	map <- genmapB36$genmap
	if (is.null(centromere)) centromere <- genmapB36$centromere else {
	  if (!is.data.frame(centromere)) stop("'centromere' argument must be NULL or a data frame")
	  if (names(centromere)!=c("chr", "start", "end")) stop("'centromere' must contains 'chr', 'start', 'end' variables, one line by chromosom")
	}
  } else if (map=='B37') {
    data('genmapB37')
	map <- genmapB37$genmap
	if (is.null(centromere)) centromere <- genmapB37$centromere else {
	  if (!is.data.frame(centromere)) stop("'centromere' argument must be NULL or a data frame")
	  if (names(centromere)!=c("chr", "start", "end")) stop("'centromere' must contains 'chr', 'start', 'end' variables, one line by chromosom")
	}
  } else if (is.data.frame(map)) {
    if (names(map)!=c("chr", "base", "rate.cM.Mb", "cM")) stop("'map' must contains 'chr', 'base', 'rate.cM.Mb' and 'cM' variables")
	if (is.null(centromere)) stop("If 'map' is a data frame, 'centromere' argument must be given")
	if (!is.data.frame(centromere)) stop("'centromere' argument must be a data frame")
	if (names(centromere)!=c("chr", "start", "end")) stop("'centromere' must contains 'chr', 'start', 'end' variables, one line by chromosom")
  } else stop("'map' must be equal to 'B36' or 'B37' or a data frame")
  
  if (is.character(bed))
  {
    snps <- read.table(paste(bed, '.bim', sep=''))
  } else if (class(bed) == "bed.matrix")
  {
    snps <- bed@snps	
  } else {
    stopCluster(cl)
	stop('"bed" must be a string or a bed.matrix.') 
  }
  temp <- order(snps$chr, snps$pos)
  if (sum(temp!=1:nrow(snps))>0)
  {
    warning('SNPs are not ordered')
	snps <- snps[temp,]
	if (is.character(bed)) bed <- read.bed.matrix(bed)[,temp] else bed <- bed[,temp]
  }
  snps$index <- 1:nrow(snps)
 
  result <- NULL
  for (i in unique(snps$chr))
  {
    w <- snps[snps$chr==i,]
    if (unit=='indices')
	{
	  temp <- data.frame(chr=i, start=NA, end=NA, index.start= w$index[seq(1, floor((nrow(w)-1)/sliding)*sliding, by=sliding)], index.end=NA)
	  temp$index.end <- sapply(temp$index.start, function(x) min(x + windows, w$index[nrow(w)]) )
	  temp$start <- snps$pos[temp$index.start]
	  temp$end <- snps$pos[temp$index.end]
	} else if (unit=='bases')
	{
	  temp <- data.frame(chr=i, start=seq(0, as.integer(ceiling( max(w$pos)/sliding )*sliding), by=sliding), end=NA, index.start=NA, index.end=NA)
	  temp$end <- temp$start + windows	 
      tt <- apply(cbind(temp$start,temp$end), 1, function(x) {
	                      ww <- ifelse( w$pos>=x[1] & w$pos<x[2], w$pos, NA )
						  if (sum(!is.na(ww))>0) return( c(which.min(ww), which.max(ww)) )
						  else return(rep(NA,2)) } )  
	  temp$index.start <- w$index[tt[1,]]
	  temp$index.end <- w$index[tt[2,]]
	} else stop("'unit' must be equal to 'base' or 'indices'")
	
	temp$recombi_snps <- apply(cbind(temp$start, temp$end), 1, function(x) mean( map$rate.cM.Mb[ map$chr==i & map$base>=x[1] & map$base<=x[2] ] ) )
	temp$recombi_region <- apply(cbind(temp$start, temp$end), 1, function(x) {
        w <- which(map$chr==i & map$base>=x[1] & map$base<=x[2]); 
		if (length(w)==0) return(NA)
		else if (length(w)==1) return(map$rate.cM.Mb[w])
		else return( diff(range( map$cM[ w ] ))/diff(range( map$base[ w ]/1e6 )) ) })

    result <- rbind(result, temp)
  }
  rm(temp, w)
  result$centro <- ( (result$start > centromere$start[result$chr] & result$start < centromere$end[result$chr]) |
                     (result$end > centromere$start[result$chr] & result$end < centromere$end[result$chr]) |
					 (result$end > centromere$end[result$chr] & result$start < centromere$start[result$chr]) )
 
  #if(!file.exists(paste(bed, 'bed',sep='.')))
  
  result$num <- NA
  result$LD <- NA
  result$LD_sd <- NA
  result$maf <- NA
  result$maf_sd <- NA
  result$H <- NA
  result$H_sd <- NA
  result$Hz <- NA
  result$Hz_sd <- NA
  
  w <- names(relatedness)
  ww <- names(FUN)
  for (i in w) for (j in ww) result[paste(i, j, sep='_')] <- NA

  if (thread>1)
  {
    if (class(bed) == "bed.matrix") {
	  time <- paste(c('temp_gaston.pop', strsplit(date(), ' ')[[1]][strsplit(date(), ' ')[[1]]!='']), collapse='_')
	  write.bed.matrix(bed, time, rds=NULL)
	  bed <- time
    }
	
    RNGkind("L'Ecuyer-CMRG")
	cl <- makeForkCluster(nnodes=thread)
    #cl <- makeCluster(thread)
	#lib <- eval(.libPaths(), envir=.GlobalEnv)
	#lib <- lib[length(lib)]
    #clusterExport(cl, c("result", "relatedness", "FUN", "LD", "bed", "couple", "lib"), envir=environment())
	#clusterCall(cl, function() .libPaths(lib))
    #clusterCall(cl, function() library(gaston))
    #clusterCall(cl, function() library(gaston.pop))
  
	result[,-(1:8)] <- matrix( parRapply(cl, cbind(result$chr, result$start, result$end, result$index.start, result$index.end), function(xx) {
	  if (is.na(xx[4]) | is.na(xx[5])) return(rep(NA, 9+length(relatedness)*length(FUN)))
	  if (is.character(bed)) x <- gaston.pop:::read.bed.matrix.part(bed, beg=xx[4], end=xx[5])
      else if (class(bed) == "bed.matrix") x <- select.snps(bed, chr==xx[1] & pos>=xx[2] & pos<xx[3])	

	  if(ncol(x)==0) return(rep(NA, 9+length(relatedness)*length(FUN)))
      if (ncol(x)>0) {
	    standardize(x) <- "mu_sigma"
        if (!is.null(LD.thin) & ncol(x)>1) x <- LD.thin(x, LD.thin, extract=TRUE)
        r <- ncol(x)
        #if (sum(map[[paste('chr', unique(x@snps$chr), sep='')]]$base %in% x@snps$pos)>0)  t$recombi <- map[[paste('chr', unique(x@snps$chr), sep='')]]( map[[paste('chr', unique(x@snps$chr), sep='')]]$rate.cM.Mb[ map[[paste('chr', unique(x@snps$chr), sep='')]]$base %in% x@snps$pos ], na.rm=TRUE)
      }

      # LD
      L <- LD(x, lim=c(1,ncol(x)))
	  L <- L[lower.tri(L, diag=FALSE)]

      r <- c(r, mean(L) )
      r <- c(r, sd(L) )
      r <- c(r, mean(x@snps$maf) )
      r <- c(r, sd(x@snps$maf) )
      r <- c(r, mean(2*x@snps$maf*(1-x@snps$maf)) )
      r <- c(r, sd(2*x@snps$maf*(1-x@snps$maf)) )
      r <- c(r, mean(x@snps$hz) )
      r <- c(r, sd(x@snps$hz) )

	  standardize(x) <- 'p'
      for (i in 1:length(relatedness))
      {
        app <- relatedness[[i]](x)
		for (j in 1:length(FUN)) {
		  t <- FUN[[j]](app)
		  if (length(t)>1) warning("FUN element(s) must be function with result value of length 1")
		  else  r <- c(r, t)
		}
      }
      return(r) }), ncol=9+length(relatedness)*length(FUN), byrow=TRUE)
    stopCluster(cl)
	file.remove(paste(time, 'bed', sep='.'))
	file.remove(paste(time, 'bim', sep='.'))
	file.remove(paste(time, 'fam', sep='.'))
  } else {
    for (k in 1:nrow(result))
	{
	  if (is.character(bed)) x <- gaston.pop:::read.bed.matrix.part(bed, beg=result$index.start[k], end=result$index.end[k])
      else if (class(bed) == "bed.matrix") x <- select.snps(bed, chr==result$chr[k] & pos>=result$start[k] & pos<result$end[k])	

	  if(ncol(x)==0) next
      if (ncol(x)>0) {
        standardize(x) <- 'mu_sigma'
        if (!is.null(LD.thin) & ncol(x)>1) x <- LD.thin(x, LD.thin, extract=TRUE)
        result$num[k] <- ncol(x)
        #if (sum(map[[paste('chr', unique(x@snps$chr), sep='')]]$base %in% x@snps$pos)>0)  t$recombi <- map[[paste('chr', unique(x@snps$chr), sep='')]]( map[[paste('chr', unique(x@snps$chr), sep='')]]$rate.cM.Mb[ map[[paste('chr', unique(x@snps$chr), sep='')]]$base %in% x@snps$pos ], na.rm=TRUE)
      }

      # LD
      L <- LD(x, lim=c(1,ncol(x)))
	  L <- L[lower.tri(L, diag=FALSE)]
      result$LD[k] <- mean(L)
      result$LD_sd[k] <- sd(L)
      result$maf[k] <- mean(x@snps$maf)
      result$maf_sd[k] <- sd(x@snps$maf)
      result$H[k] <- mean(mean(2*x@snps$maf*(1-x@snps$maf)))
      result$H_sd[k] <- sd(2*x@snps$maf*(1-x@snps$maf))
      result$Hz[k] <- mean(x@snps$hz)
      result$Hz_sd[k] <- sd(x@snps$hz)

	  standardize(x) <- 'p'
      for (i in 1:length(relatedness))
      {
        app <- relatedness[[i]](x)
		for (j in 1:length(FUN)) {
		  t <- FUN[[j]](app)
		  if (length(t)>1) warning("FUN element(s) must be function with result value of length 1")
		  else  result[k, paste(w[i], ww[j], sep='_')] <- t
		}
      }
    }
  }
  
  return(result)
}
