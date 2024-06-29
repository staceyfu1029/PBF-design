library(sm)
K <- 4
phi <- 0.25
cs <- 3

working <- dt[dt$K==K & dt$phi==phi & dt$cohortsize==cs & dt$method=="boin",]
load(paste0("BOIN","-",phi*100,"-",K,"-",cs,".RData"))

mtd <- rep(0,n.scene)
tol.mtd <- matrix(nrow = n.scene, ncol = K)
for(i in 1:n.scene){
  mtd[i] <- unname(which.min(abs(skeleton_list[i,]-phi)))
  for(j in 1:K){
    tol.mtd[i,j] <- ifelse(skeleton_list[i,j]>=phi-0.05 & skeleton_list[i,j]<=phi+0.05,1,0)
  }
}

# sm.density.compare(-working$y,mtd,xlab=paste0("Diff in PCS, K=",K," N=",n," phi=",phi))
# colfill<-2:(2+K)
# legend("topright",legend=1:K, fill=colfill)
# 
# hist(-working$y[mtd==1],breaks=30)

p <- data.frame(y=working$y,mtd=mtd)
dens <- tapply(-working$y, mtd, density)

ylim <- 0
for(i in 1:K){
  ylim <- max(ylim,dens[[i]]$y)
}
plot(density(-working$y),xlab=paste0("Diff in PCS, K=",K,", N=",n,", \u03a6=",phi),
     ylim=c(0,ceiling(ylim/10)*10),main = "",lwd=2)
for(i in 1:K){
  lines(dens[[i]],col=(i+1))
}

# Legend
legend("topright", legend = c("overall",1:K),
       lty = 1, col = 1:(K+1))

par(mfrow=c(2,2))
for(k in 1:length(K.list)){
  K <- K.list[k]
  working <- dt[dt$K==K & dt$phi==phi & dt$cohortsize==cs & dt$method=="boin",]
  load(paste0("BOIN","-",phi*100,"-",K,"-",cs,".RData"))
  
  mtd <- rep(0,n.scene)
  tol.mtd <- matrix(nrow = n.scene, ncol = K)
  for(i in 1:n.scene){
    mtd[i] <- unname(which.min(abs(skeleton_list[i,]-phi)))
    for(j in 1:K){
      tol.mtd[i,j] <- ifelse(skeleton_list[i,j]>=phi-0.05 & skeleton_list[i,j]<=phi+0.05,1,0)
    }
  }
  
  # sm.density.compare(-working$y,mtd,xlab=paste0("Diff in PCS, K=",K," N=",n," phi=",phi))
  # colfill<-2:(2+K)
  # legend("topright",legend=1:K, fill=colfill)
  # 
  # hist(-working$y[mtd==1],breaks=30)
  
  p <- data.frame(y=working$y,mtd=mtd)
  dens <- tapply(-working$y, mtd, density)
  
  ylim <- 0
  for(i in 1:K){
    ylim <- max(ylim,dens[[i]]$y)
  }
  plot(density(-working$y),xlab=paste0("Diff in PCS of PoP vs BOIN, K=",K,", N=",n,", \u03a6=",phi),
       ylim=c(0,ceiling(ylim/10)*10),main = "",lwd=2)
  for(i in 1:K){
    lines(dens[[i]],col=(i+1))
  }
  lines(density(-working$y),lwd=2)
  
  # Legend
  legend("topright", legend = c("overall",1:K),
         lty = 1, col = 1:(K+1))
}



