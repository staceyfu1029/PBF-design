rm(list = ls())

library(dfcrm)
library(BOIN)
library(Keyboard)

#### Functions ------------------------------------------

trial.pop <- function(lower,upper,skeleton,start)
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
    
    # Starting dose
    start.dose <- start
    
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
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
          break
        }
      }
      if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
        dose.elim[dose.next:K] <- 0  ## Exclude for being overly toxic
        if(sum(dose.elim==1)==0){
          early[count] <- 1
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

trial.boin <- function(lower,upper,skeleton,start)
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
    # Starting dose
    start.dose <- start
    
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
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
      if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
        dose.elim[dose.next:K] <- 0
        if(dose.elim[1]==0){
          early[count] <- 1
          #mtd[count] <- 1
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

trial.kb <- function(lower,upper,skeleton,start)
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
    # Starting dose
    start.dose <- start
    
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
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
      if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
        dose.elim[dose.next:K] <- 0
        if(dose.elim[1]==0){
          early[count] <- 1
          #mtd[count] <- 1
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

trial.crm <- function(lower,upper,skeleton,start)
{
  true.mtd <- which.min(abs(skeleton-phi))
  K <- length(skeleton)
  mtd <- rep(NA,n.sim)
  num.p <- num.tox <- matrix(nrow = n.sim,ncol = K)
  risk.over <- rep(NA,n.sim)
  risk.under <- rep(NA,n.sim)
  early <- rep(0,n.sim)
  
  prior <- getprior(halfwidth = 0.05, target = phi, nu=ceiling(K/2),nlevel = K, model = "logistic")
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
  
    level <- numeric()
    y <- numeric()
    
    while(s>0){
      dose.treated[dose.next] <- dose.treated[dose.next]+1
      level <- c(level, dose.next)
      dlt <- toxic[s,dose.next]
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      y <- c(y, dlt)
      s <- s-1
      
      if(dlt==1){
        mtd.hat <- crm(prior, phi, y, level)$mtd
        if(mtd.hat<dose.next){ # if observed dlt <= boundary, escalate
          dose.next <- max(1,dose.next-1)
        }else if(mtd.hat>dose.next){
          dose.next <- min(K,dose.next+1)
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
      level <- c(level, rep(dose.next,s.tr))
      dlt <- sum(toxic[s+1-(1:s.tr),dose.next])
      y <- c(y, toxic[s+1-(1:s.tr),dose.next])
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      s <- s-s.tr
      
      ## No exclusion rules for CRM, but an early exclusion rule is applied
      if (1 - pbeta(phi, dose.dlt[1] + 1, dose.treated[1] - dose.dlt[1] + 1) > 0.95) {
        mtd[count] <- 99
        early[count] <- 1
        break
      }
      
      ## Transition decision
      mtd.hat <- crm(prior, phi, y, level)$mtd
      
      if(mtd.hat<dose.next){ 
        dose.next <- max(1,dose.next-1)
      }else if(mtd.hat>dose.next){
        dose.next <- min(K,dose.next+1)
      }
    }
    
    
    if(is.na(mtd[count])){
      mtd[count] <- crm(prior, phi, y, level)$mtd
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

prbf01 <- function(n,y,phi){
  p <- beta(y+2,n-y+1)/beta(1+y,n-y+1)
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

#### Overall setup --------------------------------------
seed <- 2024
K <- 6
n <-  36
cohortsize <- 1
phi <- 0.25
start <- 1
n.cohort <- n/cohortsize
n.sim <- 20000
risk.cutoff <- 0.7
digits <- 1

#### run multiple scenarios ----------------------------

skeleton.list <- readxl::read_excel("prespecified_scenes_boin.xlsx",col_names = FALSE)

for(index in 1:dim(skeleton.list)[1]){
  skeleton <- round(as.numeric(skeleton.list[index,]),2)
  
  ## run pop design
  pbfbound <- bound(n,1,phi)
  lower <- pbfbound[1,]
  upper <- pbfbound[2,]
  elim.lower <- pbfbound[3,]
  elim.upper <- pbfbound[4,]
  
  set.seed(seed)
  true.mtd <- which.min(abs(skeleton-phi))
  t <- trial.pop(lower,upper,skeleton,start = start)
  
  tab <- data.frame(matrix(nrow = 3,ncol = K+4))
  tab[1,1:K] <- skeleton
  tab[2,1:(K+1)] <- round(cont(t$mtd,K)*100,digits = digits)
  tab[3,1:K] <- round(colMeans(t$num.p,na.rm = T),digits = digits) #/n*100
  tab[2,K+3] <- round(mean(t$early)*100,digits = digits)
  tab[3,K+2] <- round(mean(rowSums(t$num.p)),digits = digits)
  tab[2,K+4] <- round(mean(t$risk.over*100,na.rm = T),digits = digits)
  colnames(tab) <- c(seq(1,K),"99","Number of pts","Early stopping","Risk of overdosing")
  rownames(tab) <- c("True DLT rates","PCS","PCA")
  
  if(index==1){
    out <- cbind(design=c("","PoP",""),metrics=rownames(tab),tab)
  }else{
    out <- rbind(out,cbind(design=c("","PoP",""),metrics=rownames(tab),tab))
  }
  
  ## run BOIN design
  boin.bound <- get.boundary(target = phi,ncohort = n,cohortsize = 1)
  lower <- boin.bound$boundary_tab[2,]
  upper <- boin.bound$boundary_tab[3,]
  elim.upper <- boin.bound$boundary_tab[4,]
  elim.upper[which(is.na(elim.upper))] <- Inf
  
  set.seed(seed)
  true.mtd <- which.min(abs(skeleton-phi))
  t <- trial.boin(lower,upper,skeleton,start = start)
  
  tab <- data.frame(matrix(nrow = 3,ncol = K+4))
  tab[1,1:K] <- skeleton
  tab[2,1:(K+1)] <- round(cont(t$mtd,K)*100,digits = digits)
  tab[3,1:K] <- round(colMeans(t$num.p,na.rm = T),digits = digits) #/n*100
  tab[2,K+3] <- round(mean(t$early)*100,digits = digits)
  tab[3,K+2] <- round(mean(rowSums(t$num.p)),digits = digits)
  tab[2,K+4] <- round(mean(t$risk.over*100,na.rm = T),digits = digits)
  colnames(tab) <- c(seq(1,K),"99","Number of pts","Early stopping","Risk of overdosing")
  rownames(tab) <- c("True DLT rates","PCS","PCA")
  
  out <- rbind(out,cbind(design=c("BOIN",""),metrics=rownames(tab[-1,]),tab[-1,]))
  
  ## run Keyboard design
  kb.bound <- get.boundary.kb(target = phi,ncohort = n,cohortsize = 1)
  lower <- kb.bound$boundary_tab[2,]
  upper <- kb.bound$boundary_tab[3,]
  elim.upper <- kb.bound$boundary_tab[4,]
  elim.upper[which(is.na(elim.upper))] <- Inf
  
  set.seed(seed)
  true.mtd <- which.min(abs(skeleton-phi))
  t <- trial.boin(lower,upper,skeleton,start = start)
  
  tab <- data.frame(matrix(nrow = 3,ncol = K+4))
  tab[1,1:K] <- skeleton
  tab[2,1:(K+1)] <- round(cont(t$mtd,K)*100,digits = digits)
  tab[3,1:K] <- round(colMeans(t$num.p,na.rm = T),digits = digits) #/n*100
  tab[2,K+3] <- round(mean(t$early)*100,digits = digits)
  tab[3,K+2] <- round(mean(rowSums(t$num.p)),digits = digits)
  tab[2,K+4] <- round(mean(t$risk.over*100,na.rm = T),digits = digits)
  colnames(tab) <- c(seq(1,K),"99","Number of pts","Early stopping","Risk of overdosing")
  rownames(tab) <- c("True DLT rates","PCS","PCA")
  
  out <- rbind(out,cbind(design=c("Keyboard",""),metrics=rownames(tab[-1,]),tab[-1,]))
  
  ## run CRM design
  set.seed(seed)
  true.mtd <- which.min(abs(skeleton-phi))
  t <- trial.crm(lower,upper,skeleton,start = start)

  tab <- data.frame(matrix(nrow = 3,ncol = K+4))
  tab[1,1:K] <- skeleton
  tab[2,1:(K+1)] <- round(cont(t$mtd,K)*100,digits = digits)
  tab[3,1:K] <- round(colMeans(t$num.p,na.rm = T),digits = digits) #/n*100
  tab[2,K+3] <- round(mean(t$early)*100,digits = digits)
  tab[3,K+2] <- round(mean(rowSums(t$num.p)),digits = digits)
  tab[2,K+4] <- round(mean(t$risk.over*100,na.rm = T),digits = digits)
  colnames(tab) <- c(seq(1,K),"99","Number of pts","Early stopping","Risk of overdosing")
  rownames(tab) <- c("True DLT rates","PCS","PCA")
  out <- rbind(out,cbind(design=c("CRM",""),metrics=rownames(tab[-1,]),tab[-1,]))
}

writexl::write_xlsx(out,path = "Table 3.xlsx")
