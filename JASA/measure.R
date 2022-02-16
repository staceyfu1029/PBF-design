measure <- function(design.mtd,design.allo,design){
  pcs <- rep(0,n.scene)
  pcs.range <- rep(0,n.scene)
  pca <- rep(0,n.scene)
  pca.range <- rep(0,n.scene)
  overdose <- rep(0,n.scene)
  underdose <- rep(0,n.scene)
  pos <- rep(0,n.scene)
  for(i in 1:n.scene){
    pcs[i] <- design.mtd[i,mtd[i]]
    if(sum(tol.mtd[i,])==0){
      pcs.range[i] <- NA
      pca.range[i] <- NA
    }else{
      pcs.range[i] <- sum(design.mtd[i,which(tol.mtd[i,]==1)])
      pca.range[i] <- sum(design.allo[i,which(tol.mtd[i,]==1)])
    }
    pca[i] <- design.allo[i,mtd[i]]
    if(mtd[i]==K){
      overdose[i] <- 0
      pos[i] <- 0
    }else{
      overdose[i] <- sum(design.allo[i,(mtd[i]+1):K])
      pos[i] <- sum(design.mtd[i,(mtd[i]+1):K])
    }
    if(mtd[i]==1){
      underdose[i] <- 0
    }else{
      underdose[i] <- sum(design.allo[i,1:mtd[i]])
    }
  }
  act.ss <- rowSums(design.allo)*n
  risk.over <- design[,3]
  risk.under <- design[,4]
  d <- data.frame(pcs=pcs,pcs.range=pcs.range,pca=pca,pca.range=pca.range,
                  overdose=overdose,risk.over=risk.over,
                  underdose=underdose,risk.under=risk.under,
                  pos=pos,
                  act.ss=act.ss)
  d
}
