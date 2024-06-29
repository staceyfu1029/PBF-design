rm(list=ls())
## The following code are used to run trials for random scenarios

## load skeleton list
library(readxl)
skeleton_list <- read_excel("target15-4.xlsx", col_names = FALSE)

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
    for(i in 1:n){
      toxic[i,1] <- rbinom(1,1,prob = skeleton[1])
      for(j in 2:K){
        if(toxic[i,j-1]==1){
          toxic[i,j] <- 1
        }else{
          phi1 <- skeleton[j-1]
          phi2 <- skeleton[j]
          toxic[i,j] <- rbinom(1,1,prob = (phi2-phi1)/(1-phi1))
        }
      }
    }
    
    # }
    # Starting dose
    start.dose <- start
    
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
    dose.elim <- rep(1,K)
    
    s <- n
    
    # while(s>0){
    #   dose.treated[dose.next] <- dose.treated[dose.next]+1
    #   dlt <- toxic[s,dose.next]
    #   dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
    #   s <- s-1
    # 
    #   if(dlt==1){
    #     if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){ # if observed dlt <= boundary, escalate
    #       dose.next <- min(K,dose.next+1)
    #     }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
    #       dose.next <- max(dose.next-1,1)
    #     }
    #     break
    #   }else{
    #     dose.next <- min(K,dose.next+1)
    #   }
    # }
    
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
        if(dose.next==K){
          break
        }else{
          dose.next <- min(K,dose.next+1)
        }
      }
    }
    
    while(s>0){
      s.tr <- min(s,cohortsize)
      dose.treated[dose.next] <- dose.treated[dose.next]+s.tr
      dlt <- sum(toxic[s+1-(1:s.tr),dose.next])
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      s <- s-s.tr
      
      ## Exclusion decision
      if(dose.dlt[dose.next]<=elim.lower[dose.treated[dose.next]]){
        dose.elim[1:dose.next] <- -1 ## Exclude for being subtherapeutic
        if(sum(dose.elim==1)==0){
          early[count] <- 1
          #mtd[count] <- dose.next
          break
        }
      }
      if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
        dose.elim[dose.next:K] <- 0  ## Exclude for being overly toxic
        if(sum(dose.elim==1)==0){
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
      #mtd[count] <- select.mtd(target = phi,npts = dose.treated,ntox = dose.dlt)$MTD
      mtd[count] <- iso.pop(dose.dlt,dose.treated)
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

iso.pop <- function(p1,p0){
  l <- which(p0>0)
  p <- (p1[l]+0.05)/(p0[l]+0.1)
  p.var = (p1[l] + 0.05) * (p0[l] - p1[l] + 0.05)/((p0[l] + 
                                                      0.1)^2 * (p0[l] + 0.1 + 1))
  p.iso <- pava(p, wt = 1/p.var)
  p.iso = p.iso + (1:length(p.iso)) * 1e-10
  
  ## eliminate dose based on posterior probability
  elimi = rep(0, K)
  for (i in 1:length(l)) {
    if (1 - pbeta(phi, p1[l[i]] + 1, p0[l[i]] - p1[l[i]] + 1) > 0.95) {
      elimi[l[i]:K] = 1
      break
    }
  }
  m <- which(elimi!=1)
  
  l <- l[m]
  p.iso <- p.iso[m]
  if(length(l)==0) {return(99)}
  l[sort(abs(p.iso - phi), index.return = T)$ix[1]]
}

pava <- function(x, wt = rep(1, length(x))) {
  n <- length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}

cont <- function(x,n.level){
  ret <- rep(0,n.level+1)
  for(i in c(1:n.level)){
    ret[i] <- mean(x==i,na.rm = T)
  }
  ret[n.level+1] <- mean(x==99,na.rm = T)
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
prbf01 <- function(n,y,phi){
  #p <- beta(y+2,n-y+1)/beta(1+y,n-y+1)
  p <- exp(lbeta(y+2,n-y+1)-lbeta(1+y,n-y+1))
  dbinom(y,n,prob=phi)*exp(1)/dbinom(y,n,prob=p)
}
bound <- function(n.cohort,cohortsize,phi){
  lower <- upper <- 0
  lower.ex <- upper.ex <- 0
  for(i in 1:n.cohort){
    t <- i*cohortsize
    x <- 0:(i*cohortsize)
    y <- lapply(x,prbf01,n=(i*cohortsize),phi=phi)
    a <- 0
    for(j in 1:length(x)){
      if(y[[j]]<2.5){ 
        if(x[j]/(i*cohortsize)<phi){
          a[j] <- 1
        }else{
          a[j] <- -1
        }
      }else{
        a[j] <- 0
      }
    }
    lower[i] <- max(which(a==1))-1
    upper[i] <- min(which(a==-1))-1
    ex <- 0
    for(j in 1:length(x)){
      if(y[[j]]<5/24){
        if(x[j]/(i*cohortsize)<phi){
          ex[j] <- 1 # exclude for being subtherapeutic
        }else{
          ex[j] <- -1 # exclude for being too toxic
        }
      }else{
        ex[j] <- 0
      }
    }
    lower.ex[i] <- max(which(ex==1))-1
    upper.ex[i] <- min(which(ex==-1))-1
  }
  b <- rbind(lower,upper,lower.ex,upper.ex)
}

pbfbound <- bound(n,1,phi)
lower <- pbfbound[1,]
upper <- pbfbound[2,]
elim.lower <- pbfbound[3,]
elim.upper <- pbfbound[4,]


## run it!
set.seed(2024)
start <- 1
n.scene <- dim(skeleton_list)[1]
design.mtd <- data.frame(matrix(nrow=n.scene,ncol = K+1))
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

#save.image("PoP-15-4-1.RData")
