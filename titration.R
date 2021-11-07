## run trials with titration steps
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
    # Starting dose
    start.dose <- start
    
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- 1
    dose.elim <- rep(1,K)
    
    s <- n
    
    while(s>0){
      dose.treated[dose.next] <- dose.treated[dose.next]+1
      dlt <- rbinom(1,1,prob = skeleton[dose.next])
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
      dlt <- rbinom(1,s.tr,prob = skeleton[dose.next])
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
