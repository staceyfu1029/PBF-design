###



## Simulate the true toxicity profiles
##   Randomly generate scenarios by the pseudouniform algorithm 
##        (see Clertant and O'Quigley,2017)
scene <- function(phi,K){ 
  MTD <- sample(K,1)
  M <- rbeta(1,max(K-MTD,0.5),1)
  B <- phi+(1-phi)*M
  skeleton <- rep(0,K)
  if(MTD == 1){
    skeleton[2:K] <- sort(runif(K-1,phi,B))
    skeleton[1] <- runif(1,max(0,2*phi-skeleton[2]),skeleton[2])
  }else if(MTD == K){
    skeleton[1:(K-1)] <- sort(runif(K-1,0,phi))
    skeleton[K] <- runif(1,skeleton[K-1],min(2*phi-skeleton[K-1],B))
  }else{
    skeleton[1:(MTD-1)] <- sort(runif(MTD-1,0,phi))
    skeleton[(MTD+1):K] <- sort(runif(K-MTD,phi,B))
    d <- min(phi-skeleton[MTD-1],skeleton[MTD+1]-phi)
    skeleton[MTD] <- runif(1,phi-d,phi+d)
  }
  round(skeleton,digits = 2)
}


## run simulations
## run trials including EARLY TERMINATION rules.
trial <- function(lower,upper,skeleton,start)
{
  true.mtd <- which.min(abs(skeleton-phi))
  K <- length(skeleton)
  mtd <- rep(NA,n.sim)
  num.p <- matrix(nrow = n.sim,ncol = K)
  risk.over <- rep(NA,n.sim)
  risk.under <- rep(NA,n.sim)
  early <- rep(0,n.sim)
  for(count in 1:n.sim){
    if(start==1){
      start.dose <- 1
    }else if(start==2){
      start.dose <- sample(c(ceiling(K/2),ceiling(K/2+0.5)),1)
    }else{
      start.dose <- sample(1:K,1)
    }
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
    dose.elim <- rep(1,K)
    
    for(i in 1:n.cohort){
      dose.treated[dose.next] <- dose.treated[dose.next]+1
      dlt <- rbinom(1,cohortsize,prob = skeleton[dose.next])
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      if(dose.dlt[dose.next]<=elim.lower[dose.treated[dose.next]]){
        dose.elim[1:dose.next] <- 0
        if(sum(dose.elim)==0){
          early[count] <- 1
          mtd[count] <- dose.next
          break
        }
      }
      if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
        dose.elim[dose.next:K] <- 0
        if(sum(dose.elim)==0){
          early[count] <- 1
          mtd[count] <- dose.next
          break
        }
      }
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
      mtd[count] <- iso(dose.dlt,dose.treated*cohortsize)
    }
    if(is.na(mtd[count])){
      next
    }else{
      num.p[count,] <- dose.treated*cohortsize
      risk.over[count] <- 0
      risk.under[count] <- 0
      if(true.mtd==1){
        if(sum(dose.treated[2:K])*cohortsize>risk.cutoff*n){
          risk.over[count] <- 1
        }
      }else if(true.mtd==K){
        if(sum(dose.treated[1:(K-1)])*cohortsize>risk.cutoff*n){
          risk.under[count] <- 1
        }
      }else{
        if(sum(dose.treated[1:(true.mtd-1)])*cohortsize>risk.cutoff*n){
          risk.under[count] <- 1
        }else if(sum(dose.treated[(true.mtd+1):K])*cohortsize>risk.cutoff*n){
          risk.over[count] <- 1
        }
      }
    }
  }
  return(list(num.p=num.p,mtd=mtd,early=early,
              risk.over=risk.over,risk.under=risk.under))
}




## Simulate trials, without early termination rules
trial <- function(lower,upper,skeleton)
{
  true.mtd <- which.min(abs(skeleton-phi))
  K <- length(skeleton)
  start.dose <- 1
  mtd <- rep(0,n.sim)
  num.p <- rep(0,K)
  risk.over <- rep(0,n.sim)
  risk.under <- rep(0,n.sim)
  for(count in 1:n.sim){
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
    
    for(i in 1:n.cohort){
      dose.treated[dose.next] <- dose.treated[dose.next]+1
      dlt <- rbinom(1,cohortsize,prob = skeleton[dose.next])
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){
        dose.next <- dose.next+1
      }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
        dose.next <- dose.next-1
      }
      if(dose.next<1){
        dose.next <- 1
      }else if(dose.next>length(skeleton)){
        dose.next <- length(skeleton)
      }
    }
    mtd[count] <- iso(dose.dlt,dose.treated*cohortsize)
    num.p <- num.p+dose.treated*cohortsize
    if(true.mtd==1){
      if(sum(dose.treated[2:K])*cohortsize>risk.cutoff*n){
        risk.over[count] <- 1
      }
    }else if(true.mtd==K){
      if(sum(dose.treated[1:(K-1)])*cohortsize>risk.cutoff*n){
        risk.under[count] <- 1
      }
    }else{
      if(sum(dose.treated[1:(true.mtd-1)])*cohortsize>risk.cutoff*n){
        risk.under[count] <- 1
      }else if(sum(dose.treated[(true.mtd+1):K])*cohortsize>risk.cutoff*n){
        risk.over[count] <- 1
      }
    }
  }
  return(list(num.p=num.p,mtd=mtd,
              risk.over=risk.over,risk.under=risk.under))
}
