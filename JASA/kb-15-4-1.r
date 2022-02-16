rm(list=ls())
## The following code are used to run trials for random scenarios

## load skeleton list
library(readxl)
skeleton_list <- read_excel("/storage/home/cxf637/pbf/target15-4.xlsx", 
                            col_names = FALSE)

## run trials
trial <- function(lower,upper,skeleton,start)
{
  true.mtd <- which.min(abs(skeleton-phi))
  K <- length(skeleton)
  mtd <- rep(NA,n.sim)
  num.p <- num.tox <- matrix(nrow = n.sim,ncol = K)
  risk.over <- rep(NA,n.sim)
  risk.under <- rep(NA,n.sim)
  early <- rep(0,n.sim)
  
  for(count in 1:n.sim){
    toxic <- matrix(nrow = n,ncol = K)
    for(j in 1:K){
      toxic[,j] <- rbinom(n,1,prob = skeleton[j])
    }
    # Starting dose
    start.dose <- start
    
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- 1
    dose.elim <- rep(1,K)
    
    s <- n
    
    while(s>0){
      dose.treated[dose.next] <- dose.treated[dose.next]+1
      dlt <- toxic[s,dose.next]
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      s <- s-1
      
      if(dlt==1){
        if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){ # if observed dlt <= boundary, escalate
          dose.next <- min(K,dose.next+1)
        }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
          dose.next <- max(dose.next-1,1)
        }
        break
      }else{
        dose.next <- min(K,dose.next+1)
      }
    }
    
    while(s>0){
      s.tr <- min(s,cohortsize)
      dose.treated[dose.next] <- dose.treated[dose.next]+s.tr
      dlt <- sum(toxic[s+1-(1:s.tr),dose.next])
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      s <- s-s.tr
      
      ## Exclusion decision
      if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
        dose.elim[dose.next:K] <- 0
        if(dose.elim[1]==0){
          early[count] <- 1
          #mtd[count] <- dose.next
          break
        }
      }
      
      ## Transition decision
      if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){
        if(dose.next < K){
          if(dose.elim[dose.next+1]==1){
            dose.next <- dose.next+1
          }
        } 
      }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
        if(dose.next > 1){
          if(dose.elim[dose.next-1]==1){
            dose.next <- dose.next-1
          }
        }
      }
    }
    
    
    if(is.na(mtd[count])){
      #mtd[count] <- select.mtd(target = phi,npts = dose.treated*cohortsize,ntox = dose.dlt)$MTD
      mtd[count] <- iso(dose.dlt,dose.treated)
    }
    if(is.na(mtd[count])){
      next
    }else{
      num.p[count,] <- dose.treated
      num.tox[count,] <- dose.dlt
      risk.over[count] <- 0
      risk.under[count] <- 0
      if(true.mtd==1){
        if(sum(dose.treated[2:K])>risk.cutoff*n){
          risk.over[count] <- 1
        }
      }else if(true.mtd==K){
        if(sum(dose.treated[1:(K-1)])>risk.cutoff*n){
          risk.under[count] <- 1
        }
      }else{
        if(sum(dose.treated[1:(true.mtd-1)])>risk.cutoff*n){
          risk.under[count] <- 1
        }else if(sum(dose.treated[(true.mtd+1):K])>risk.cutoff*n){
          risk.over[count] <- 1
        }
      }
    }
  }
  return(list(num.p=num.p,mtd=mtd,early=early,num.tox=num.tox,
              risk.over=risk.over,risk.under=risk.under))
}

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

cont <- function(x,n.level){
  ret <- rep(0,n.level)
  for(i in 1:n.level){
    ret[i] <- mean(x==i,na.rm = T)
  }
  ret
}

## simulation setting
phi <- 0.15
K <- 4
n <- 30
cohortsize <- 1
n.cohort <- n/cohortsize
n.sim <- 20000
risk.cutoff <- 0.7


## BOIN boundary
library(BOIN)
library(Keyboard)

kb.bound <- get.boundary.kb(target = phi,ncohort = n,cohortsize = 1)
lower <- kb.bound$boundary_tab[2,]
upper <- kb.bound$boundary_tab[3,]
elim.upper <- kb.bound$boundary_tab[4,]
elim.upper[which(is.na(elim.upper))] <- Inf


## run it!
set.seed(1212)
start <- 1
n.scene <- dim(skeleton_list)[1]
design.mtd <- data.frame(matrix(nrow=n.scene,ncol = K))
design.allo <- data.frame(matrix(nrow=n.scene,ncol = K))
design <- matrix(nrow=n.scene,ncol = 5)
colnames(design) <- c("pcs","pca","risk.over","risk.under","earlyrate") 
Start.time <- Sys.time()
for(i in 1:n.scene){
  skeleton <- as.numeric(skeleton_list[i,])
  true.mtd <- which.min(abs(skeleton-phi))
  t <- trial(lower,upper,skeleton,start = start)
  design[i,] <- c(mean(t$mtd==true.mtd,na.rm=TRUE),mean(t$num.p[,true.mtd],na.rm = T)/n,
                  mean(t$risk.over,na.rm = T),mean(t$risk.under,na.rm = T),mean(t$early))
  design.mtd[i,] <- cont(t$mtd,K)
  design.allo[i,] <- colMeans(t$num.p,na.rm = T)/n
}
end <- Sys.time()
end-Start.time

save.image("/storage/home/cxf637/pbf/at/output/kb-15-4-1.RData")