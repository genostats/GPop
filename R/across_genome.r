app.genome <- function(bed, couple, FUN=mean, relatedness=GRM, windows, unit=c("bases", "indices"), sliding, thread=22, centromere=NULL, map='B36', LD.thin=NULL)
{
  ## Verif
  if (!is.function(FUN) & !is.list(FUN)) stop('"FUN" argument must be a function of a list of functions')
  if (is.list(FUN)) for (i in 1:length(FUN)) if (!is.function(FUN[i])) stop('"FUN" argument must be a function of a list of functions')
  if (is.function(FUN)) { f <- as.character(substitute(FUN)); FUN <- list(FUN); names(FUN) <- f; }
  
  if (!is.function(relatedness)) stop('"relatedness" argument must be a function')
  if (is.list(relatedness)) for (i in 1:length(relatedness)) if (!is.function(relatedness[i])) stop(relatedness)
  if (is.function(relatedness)) { f <- as.character(substitute(relatedness)); relatedness <- list(relatedness); names(relatedness) <- f; }

  # load of genetic map
  if (map=='B36') { data('genmapB36'); map <- genmapB36$genmap; centro <- genmapB36$centromere; }
  else if (map=='B37') { data('genmapB37'); map <- genmapB37$genmap; centro <- genmapB37$centromere; }
  else if (is.data.frame(map)) {
    if (names(map)!=c("chr", "base", "rate.cM.Mb", "cM")) stop("'map' must contains 'chr', 'base', 'rate.cM.Mb' and 'cM' variables")
	if (is.null(centromere)) stop("If 'map' is a data frame, 'centromere' argument must be given")
	if (!is.data.frame(centromere)) stop("'centromere' argument must be a data frame")
	if (names(map)!=c("chr", "start", "end")) stop("'map' must contains 'chr', 'start', 'end' variables, one line by chromosome")
  }
  else stop("'map' must be equal to 'B36' or 'B37' or a data frame")
  
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
 
  result <- NULL
  for (i in unique(snps$chr))
  {
    w <- snps[snps$chr==i,]
    if (unit=='indices')
	{
	  temp <- data.frame(chr=i, start=NA, end=NA, index.start= seq(1, floor((nrow(w)-1)/sliding)*sliding, by=sliding), index.end=NA)
	  temp$index.end <- temp$index.start + windows
	  temp$index.end[nrow(temp)] <- nrow(w)
	  temp$start <- w$pos[temp$index.start]
	  temp$end <- w$pos[temp$index.end]
	} else if (unit=='base')
	{
	  temp <- data.frame(chr=i, start=seq(0, as.integer(ceiling( max(w$pos)/sliding )*sliding), by=sliding), end=NA, index.start=NA, index.end=NA)
	  temp$end <- temp$start + windows	 
      tt <- apply(cbind(temp$start,temp$end), 1, function(x) {
	                      ww <- ifelse( w$pos>=x[1] & w$pos<x[2], w$pos, NA )
						  if (sum(!is.na(ww))>0) return( c(which.min(ww), which.max(ww)) )
						  else return(rep(NA,2)) } )  
	  temp$index.start <- tt[1,]
	  temp$index.end <- tt[2,]
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
  result$centro <- ( (result$start > centro$start[result$chr] & result$start < centro$end[result$chr]) | (result$end > centro$start[result$chr] & result$end < centro$end[result$chr]) | (result$end > centro$end[result$chr] & result$start < centro$start[result$chr]) )
 
  #if(!file.exists(paste(bed, 'bed',sep='.')))
  
  result$num <- NA
  result$LD <- NA
  result$LD_sd <- NA
  result$maf <- NA
  result$maf_sd <- NA
  
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
  
	result[,-(1:8)] <- matrix( parRapply(cl, cbind(result$chr, result$start, result$end, result$index.start, result$index.end), function(xx) {
	  if (is.na(xx[4]) | is.na(xx[5])) return(rep(NA, 5+length(relatedness)*length(FUN)))
	  if (is.character(bed)) x <- gaston.pop:::read.bed.matrix.part(bed, beg=xx[4], end=xx[5])
      else if (class(bed) == "bed.matrix") x <- select.snps(bed, chr==xx[1] & pos>=xx[2] & pos<xx[3])	

	  if(ncol(x)==0) return(rep(NA, 5+length(relatedness)*length(FUN)))
      if (ncol(x)>0) {
	    standardize(x) <- "mu_sigma"
        if (!is.null(LD.thin) & ncol(x)>1) l <- LD.thin(x, LD.thin, extract=TRUE)
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

	  standardize(x) <- 'p'
      for (i in 1:length(relatedness))
      {
        app <- relatedness[[i]](x)
		for (j in 1:length(FUN)) {
		  t <- app.pairs(app, couple, FUN=FUN[[j]])
		   if (length(t)>1) warning("FUN element(s) must be function with result value of length 1")
		  else  r <- c(r, t)
		}
      }
      return(r) }), ncol=5+length(relatedness)*length(FUN), byrow=TRUE)
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
        if (!is.null(LD.thin) & ncol(x)>1) l <- LD.thin(x, LD.thin, extract=TRUE)
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

	  standardize(x) <- 'p'
      for (i in 1:length(relatedness))
      {
        app <- relatedness[[i]](x)
		for (j in 1:length(FUN)) {
		  t <- app.pairs(app, couple, FUN=FUN[[j]])
		  if (length(t)>1) warning("FUN element(s) must be function with result value of length 1")
		  else  result[k, paste(w[i], ww[j], sep='_')] <- t
		}
      }
    }
  }
  
  return(result)
}
