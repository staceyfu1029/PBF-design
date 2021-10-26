iso <- function(p1,p0){
  l <- which(p0>0)
  p <- p1[l]/p0[l]
  if(sum(p)==0){
    return(max(l))
  }
  iso.model <- isoreg(p)
  p.iso <- fit.isoreg(iso.model,1:length(l))
  d <- abs(p.iso-phi)
  l[max(which(d==min(d)))]
}

fit.isoreg <- function(iso, x0)
{
  if(length(x0)==1){
    return(iso$yf)
  }
  o = iso$o
  if (is.null(o))
    o = 1:length(x0)
  x = unique(iso$x[o])
  y = iso$yf
  ind = cut(x0, breaks = x, labels = FALSE, include.lowest = TRUE)
  min.x <- min(x)
  max.x <- max(x)
  adjusted.knots <- iso$iKnots[c(which(iso$yf[iso$iKnots] > 0))]
  fits = sapply(seq(along = x0), function(i) {
    j = ind[i]
    
    # Find the upper and lower parts of the step
    upper.step.n <- min(which(adjusted.knots > j))
    upper.step <- adjusted.knots[upper.step.n]
    lower.step <- ifelse(upper.step.n==1, 1, adjusted.knots[upper.step.n -1] )
    
    # Perform a liner interpolation between the start and end of the step
    denom <- x[upper.step] - x[lower.step]
    denom <- ifelse(denom == 0, 1, denom)
    val <- y[lower.step] + (y[upper.step] - y[lower.step]) * (x0[i] - x[lower.step]) / (denom)
  })
  fits
}
