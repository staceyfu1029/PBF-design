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
