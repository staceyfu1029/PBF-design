## The following codes are used to generate pseudo-uniform random scenarios for 
##    simulation study.


## Input:
##    phi: target toxicity rate
##      K: dose levels
## Output:
##    one skeleton
scene <- function(phi,K){ 
  MTD <- sample(K,1)
  M <- rbeta(1,max(K-MTD,0.5),1)
  B <- phi+(1-phi)*M
  while(1){
    skeleton <- sort(runif(K,0,B))
    if(which.min(abs(skeleton-phi))==MTD){
      if(sum(skeleton >= phi-0.05 & skeleton <= phi+0.05)>0){
        break
      }
    }
  }
  skeleton
}

set.seed(1212)
#phi <- 0.40
K <- 3
n.scene <- 10000
skeleton.list <- data.frame(matrix(nrow = n.scene, ncol = K))

start <- Sys.time()
i <- 1
while(i <= n.scene){
  skeleton <- scene(phi,K)
  skeleton.list[i,] <- skeleton
  i <- i+1
}

mtd <- apply(skeleton.list,1,function(s){
  which.min(abs(s-phi))
})
table(mtd)
end <- Sys.time()
end-start

writexl::write_xlsx(skeleton.list,paste0("target",phi*100,"-",K,".xlsx"),
            col_names = FALSE)
