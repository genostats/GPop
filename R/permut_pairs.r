app.permut <- function(app, permut.pairs, group=rep(1:nrow(app)), FUN, unpairs=FALSE, fixed.pairs=NULL, B=1e5, thread=1)
{
  ## Verif
  if (!is.function(FUN)) stop('FUN must be a function')
  
  if (!is.matrix(app)) stop('app must be a matrix')
  if (ncol(app)!=nrow(app)) stop('app must be a squarred matrix')
  n <- nrow(app)

  if(!is.logical(unpairs)) stop('unpairs must be a boolean')
  
  if (!is.vector(permut.pairs) & !is.factor(permut.pairs)) stop('permut.pairs must be a vector')
  if (length(permut.pairs)!=n) stop('Dimensions of app and permut.pairs mismatch')
  
  if (!is.vector(group) & !is.factor(group)) stop('group must be a vector')
  if (length(group)!=n) stop('Dimensions of app and group mismatch')
  
  if (!is.list(fixed.pairs) & !is.null(fixed.pairs)) fixed.pairs <- list(fixed.pairs)
  
  if (is.list(fixed.pairs)) tt <- c(list(permut.pairs), fixed.pairs) else tt <- permut.pairs
  if (unpairs) app_pairs <- app_unpairs(app, tt, FUN=FUN) else app_pairs <- app_pairs(app, tt, FUN=FUN)
  
  t <- rep(NA, n)
  
  
  if (thread>1)
  {
    RNGkind("L'Ecuyer-CMRG")
    cl <- makeForkCluster(thread)
    clusterSetRNGStream(cl, runif(1)*(2**31-1) )
	
    r <- parLapply(cl, rep(floor(B/thread),thread), function(n) replicate(n, {
	           for (i in unique(group)) t[which(group==i)] <- sample(permut.pairs[which(group==i)], replace=FALSE)
			   if (is.list(fixed.pairs)) tt <- c(list(t), fixed.pairs) else tt <- t
			   if (unpairs) z <- app_unpairs(app, tt, FUN=FUN) else z <- app_pairs(app, tt, FUN=FUN)
               return(z)  } ))
	if (is.vector(r[[1]])) r <- unlist(r)
	if (is.matrix(r[[1]])) r <- do.call(cbind, r)
    stopCluster(cl)
  } else r <- replicate(B, {
               for (i in unique(group)) t[which(group==i)] <- sample(permut.pairs[which(group==i)], replace=FALSE)
			   if (is.list(fixed.pairs)) tt <- c(list(t), fixed.pairs) else tt <- t
			   if (unpairs) z <- app_unpairs(app, tt, FUN=FUN) else z <- app_pairs(app, tt, FUN=FUN)
               return(z)  } )
			   
  return(list(app_pairs=app_pairs, permut=r, twosided.p = mean(abs(r-mean(r))>abs(app_pairs-mean(r))),
              low.onesided.p=mean(r<app_pairs), up.onesided.p=mean(r>app_pairs)))
}