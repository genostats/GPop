app.pairs <- function(app, pairs, FUN=identity)
{
  ## Verif
  if (!is.function(FUN)) stop('FUN must be a function')
  
  if (!is.matrix(app)) stop('app must be a matrix')
  if (ncol(app)!=nrow(app)) stop('app must be a squarred matrix')
  n <- nrow(app)
  
   if ( (is.vector(pairs) | is.factor(pairs) ) & !is.list(pairs))
  {
    if (length(pairs)!=n) stop('Dimensions of app and pairs mismatch')
	
    temp <- tcrossprod(as.matrix(model.matrix(~-1+as.factor(pairs))))==1
    app_couple <- app[lower.tri(app, diag = FALSE) & temp ]
  } else if (is.matrix(pairs))
  {
    if (ncol(pairs)!=nrow(pairs)) stop('pairs must be a vector or a squarred matrix')
    if (nrow(app)!=nrow(pairs)) stop('Dimensions of app and pairs mismatch')
    
	app_couple <- app[lower.tri(app, diag = FALSE) & pairs ]
  } else if (is.list(pairs))
  {
    if ( !any(sapply(list(pairs, list()), function(x) (!is.vector(x) & !is.matrix(x)) | is.list(x))) ) stop('pairs must be a list of vectors or squarred matrix')
	
	temp <- matrix(TRUE, n, n)
	for (i in 1:length(pairs))
	{
	  if (is.vector(pairs[[i]]) & !is.list(pairs[[i]])) temp <- temp & tcrossprod(as.matrix(model.matrix(~-1+as.factor(pairs[[i]]))))==1
	  if (is.matrix(pairs[[i]])) temp <- temp & pairs[[i]]
	}
	app_couple <- app[lower.tri(app, diag = FALSE) & temp ]
  } else stop('pairs must be a vector, a squarred matrix or a list')
  
  return(FUN(app_couple))
}

app.unpairs <- function(app, pairs, FUN=identity)
{
  ## Verif
  if (!is.function(FUN)) stop('FUN must be a function')
  
  if (!is.matrix(app)) stop('app must be a matrix')
  if (ncol(app)!=nrow(app)) stop('app must be a squarred matrix')
  n <- nrow(app)
  
  if ( (is.vector(pairs) | is.factor(pairs) ) & !is.list(pairs))
  {
    if (length(pairs)!=n) stop('Dimensions of app and pairs mismatch')
	
    temp <- tcrossprod(as.matrix(model.matrix(~-1+as.factor(pairs))))==1
    app_couple <- app[lower.tri(app, diag = FALSE) & !temp ]
  } else if (is.matrix(pairs))
  {
    if (ncol(pairs)!=nrow(pairs)) stop('pairs must be a vector or a squarred matrix')
    if (nrow(app)!=nrow(pairs)) stop('Dimensions of app and pairs mismatch')
    
	app_couple <- app[lower.tri(app, diag = FALSE) & !pairs ]
  } else if (is.list(pairs))
  {
    if ( !any(sapply(list(pairs, list()), function(x) (!is.vector(x) & !is.matrix(x)) | is.list(x))) ) stop('pairs must be a list of vectors or squarred matrix')
	
	temp <- matrix(FALSE, n, n)
	for (i in 1:length(pairs))
	{
	  if (is.vector(pairs[[i]]) & !is.list(pairs[[i]])) temp <- temp | tcrossprod(as.matrix(model.matrix(~-1+as.factor(pairs[[i]]))))==1
	  if (is.matrix(pairs[[i]])) temp <- temp | pairs[[i]]
	}
	app_couple <- app[lower.tri(app, diag = FALSE) & !temp ]
  } else stop('pairs must be a vector, a squarred matrix or a list')
  
  return(FUN(app_couple))
}

zscore <- function(app_pairs, m, sd)
{
  ## Verif
  if (!is.vector(app_pairs)) stop('pairs must be a vector')
  
  ## relatedness between pairs of interest
  p <- mean(app_pairs, na.rm=TRUE)

  ## z-score
  z <- ( p-m )/sd
  return(z)
}