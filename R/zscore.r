zscore <- function(app_pairs, m, sd)
{
  ## Verif
  if (!is.vector(app_pairs)) stop('pairs must be a vector')
  
  ## relatedness between pairs of interest
  p <- mean(app_pairs, na.rm=TRUE))

  ## z-score
  z <- ( p-m) )/sd
  return(z)
}