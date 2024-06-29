rm(list=ls())

# Libraries

# Data
source("measure.R")

# Import R workspaces
n <- 36
K <- 6
cs <- 3
phi <- 0.20
phi.list <- c(0.2,0.25,0.3,0.15)
K.list <- c(4,6,7,9)
n.list <- c(30,36,60,96)
cs.list <- c(1,3)
method.list <- c("PoP","BOIN","Keyboard")  #"CRM"

dt <- NULL
for(phi in phi.list){
  for(k in 1:length(K.list)){
    K <- K.list[k]
    for(m in method.list){
      load(paste0(m,"-",phi*100,"-",K,"-",n.list[k],"-",cs,".RData"))
      if(m==method.list[1]){
        mtd <- rep(0,n.scene)
        for(i in 1:n.scene){
          mtd[i] <- unname(which.min(abs(skeleton_list[i,]-phi)))
        }
      }
      out <- measure(design.mtd,design.allo,design)
      if(m==method.list[1]){
        if(k==1){
          tab <- c(paste0("K=",K,",N=",n), m, format(colMeans(out[,-9],na.rm = T)*100,digits = 1,nsmall=1),mean(out[,9],na.rm = T))
        }else{
          temp <- c(paste0("K=",K,",N=",n), m, format(colMeans(out[,-9],na.rm = T)*100,digits = 1,nsmall=1),mean(out[,9],na.rm = T))
          tab <- cbind(tab,temp)
        }
      }else{
        temp <- c(NA, m, format(colMeans(out[,-9],na.rm = T)*100,digits = 1,nsmall=1),mean(out[,9],na.rm = T))
        tab <- cbind(tab,temp)
      }
      K <- K.list[k]
    }
  }
  tab <- cbind(c(paste0("phi=",phi),NA,colnames(out)),tab)
  if(is.null(dt)){
    dt <- tab
  }else{
    dt <- rbind(dt,tab)
  }
}

dt <- data.frame(dt)

writexl::write_xlsx(dt,path = "temp.xlsx",col_names = FALSE)
