## This function is used to calculate the operating characteristics of the design from the simulation outputs.

measure <- function(design.mtd,design.allo,design){
  pcs <- rep(0,n.scene)
  pca <- rep(0,n.scene)
  pos <- rep(0,n.scene) # over-dose selected
  pus <- rep(0,n.scene) # under-dose selected
  poa <- rep(0,n.scene) # over-dose treated
  for(i in 1:n.scene){
    pcs[i] <- design.mtd[i,mtd[i]]
    pca[i] <- design.allo[i,mtd[i]]
    if(mtd[i]==K){
      poa[i] <- 0
      pos[i] <- 0
    }else{
      poa[i] <- sum(design.allo[i,(mtd[i]+1):K])
      pos[i] <- sum(design.mtd[i,(mtd[i]+1):K])
    }
    if(mtd[i]==1){
      pus[i] <- 0
    }else{
      pus[i] <- sum(design.mtd[i,1:(mtd[i]-1)])
    }
  }
  pns <- design.mtd[,K+1]
  act.ss <- rowSums(design.allo)*n
  risk.over <- design[,3]
  risk.under <- design[,4]
  d <- data.frame(pcs=pcs,pca=pca,
                  pos=pos,poa=poa,risk.over=risk.over,
                  pus=pus,pns=pns,risk.under=risk.under,
                  act.ss=act.ss)
  d
}
