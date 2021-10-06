#### Packages
library(BOIN)
library(Iso)
library(writexl)
library(readxl)

#### Functions
prbf01 <- function(n,y,phi){
  p <- beta(y+2,n-y+1)/beta(1+y,n-y+1)
  dbinom(y,n,prob=phi)*exp(1)/dbinom(y,n,prob=p)
}

bound <- function(n.cohort,cohortsize,phi){
  lower <- upper <- 0
  for(i in 1:n.cohort){
    t <- i*cohortsize
    x <- 0:(i*cohortsize)
    y <- lapply(x,prbf01,n=(i*cohortsize),phi=phi)
    a <- 0
    for(j in 1:length(x)){
      if(y[[j]]<exp(1)*(1/t*log(1+t))^(1/(2*t))){ #e,exp(1)*(1/t*log(1+t))^(1/(2*t)),exp(1)*(1/t*log(1+t))^(1/(t))
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
  }
  b <- rbind(lower,upper)
}

cont <- function(x,n.level){
  ret <- rep(0,n.level)
  for(i in 1:n.level){
    ret[i] <- sum(x==i)
  }
  ret
}

iso <- function(p1,p0){
  l <- which(p0>0)
  p <- p1[l]/p0[l]
  p.iso <- pava(p)
  d <- abs(p.iso-phi)
  l[max(which(d==min(d)))]
}


#### settings

seed <- 317
#seed <- 330

skeleton_list <- as.matrix(read_excel("D:/study/Dissertation/PrBF design/Figure/Skeleton list.xlsx",
                                      col_names = FALSE))

K <- dim(skeleton_list)[2]
phi <- 0.25
n <- 30 # maximum sample size
cohortsize <- 1
n.cohort <- n/cohortsize
n.sim <- 10000
start.dose <- 4
risk.cutoff <- 0.8

#### Boundaries

pbf.bound <- bound(n.cohort,cohortsize,phi)
pbf.lower <- pbf.bound[1,]
pbf.upper <- pbf.bound[2,]

boin.bound <- get.boundary(target = phi,ncohort = n.cohort,cohortsize = cohortsize)
boin.lower <- boin.bound$boundary_tab[2,]
boin.upper <- boin.bound$boundary_tab[3,]

#### Simulation output
tab <- data.frame(matrix(nrow = 2,ncol = K+6))
tab[1,] <- c("Design",NA,"Dose level",rep(NA,5),"risk of overdosing", "risk of underdosing",
             "# Escalate","# De-escalate")
tab[2,] <- c(NA,NA,1:K,NA,NA,NA,NA)

#### Simulation
for(s in 1:dim(skeleton_list)[1]){
  skeleton <- skeleton_list[s,]
  true.mtd <- which.min(abs(skeleton-phi))
  
  ## PBF design
  set.seed(seed)
  
  mtd.pbf <- rep(0,n.sim)
  mtd.next.pbf <- rep(0,n.sim)
  num.p.pbf <- rep(0,K)
  risk.over.pbf <- rep(0,n.sim)
  risk.under.pbf <- rep(0,n.sim)
  dose.es <- rep(0,n.sim)
  dose.de <- rep(0,n.sim)
  for(count in 1:n.sim){
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
    trans <- rep(0,n.cohort-1)
    
    for(i in 1:n.cohort){
      dose.treated[dose.next] <- dose.treated[dose.next]+1
      dlt <- rbinom(1,cohortsize,prob = skeleton[dose.next])
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      if(dose.dlt[dose.next]<=pbf.lower[dose.treated[dose.next]]){
        dose.next <- dose.next+1
        trans[i] <- 1
      }else if(dose.dlt[dose.next]>=pbf.upper[dose.treated[dose.next]]){
        dose.next <- dose.next-1
        trans[i] <- -1
      }
      if(dose.next<1){
        dose.next <- 1
        trans[i] <- 0
      }else if(dose.next>length(skeleton)){
        dose.next <- length(skeleton)
        trans[i] <- 0
      }
    }
    # temp <- which.min(abs(dose.dlt/(dose.treated*cohortsize)-phi))
    iso(dose.dlt,dose.treated*cohortsize)
    mtd.pbf[count] <- iso(dose.dlt,dose.treated*cohortsize)
    mtd.next.pbf[count] <- dose.next
    dose.es[count] <- sum(trans==1)
    dose.de[count] <- sum(trans==-1)
    num.p.pbf <- num.p.pbf+dose.treated
    if(true.mtd==1){
      if(sum(dose.treated[2:K])*cohortsize>risk.cutoff*n){
        risk.over.pbf[count] <- 1
      }
    }else if(true.mtd==K){
      if(sum(dose.treated[1:(K-1)])*cohortsize>risk.cutoff*n){
        risk.under.pbf[count] <- 1
      }
    }else{
      if(sum(dose.treated[1:(true.mtd-1)])*cohortsize>risk.cutoff*n){
        risk.under.pbf[count] <- 1
      }else if(sum(dose.treated[(true.mtd+1):K])*cohortsize>risk.cutoff*n){
        risk.over.pbf[count] <- 1
      }
    }
  }
  
  ## BOIN design
  mtd.boin <- rep(0,n.sim)
  mtd.next.boin <- rep(0,n.sim)
  num.p.boin <- rep(0,K)
  risk.over.boin <- rep(0,n.sim)
  risk.under.boin <- rep(0,n.sim)
  dose.es.boin <- rep(0,n.sim)
  dose.de.boin <- rep(0,n.sim)
  set.seed(seed)
  for(count in 1:n.sim){
    dose.treated <- rep(0,K)
    dose.dlt <- rep(0,K)
    dose.next <- start.dose
    
    for(i in 1:n.cohort){
      dose.treated[dose.next] <- dose.treated[dose.next]+1
      dlt <- rbinom(1,cohortsize,prob = skeleton[dose.next])
      dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
      if(dose.dlt[dose.next]<=boin.lower[dose.treated[dose.next]]){
        dose.next <- dose.next+1
        trans[i] <- 1
      }else if(dose.dlt[dose.next]>=boin.upper[dose.treated[dose.next]]){
        dose.next <- dose.next-1
        trans[i] <- -1
      }
      if(dose.next<1){
        dose.next <- 1
        trans[i] <- 0
      }else if(dose.next>length(skeleton)){
        dose.next <- length(skeleton)
        trans[i] <- 0
      }
    }
    # temp <- which.min(abs(dose.dlt/(dose.treated*cohortsize)-phi))
    mtd.boin[count] <- iso(dose.dlt,dose.treated*cohortsize)
    mtd.next.boin[count] <- dose.next
    dose.es.boin[count] <- sum(trans==1)
    dose.de.boin[count] <- sum(trans==-1)
    num.p.boin <- num.p.boin+dose.treated
    if(true.mtd==1){
      if(sum(dose.treated[2:K])*cohortsize>risk.cutoff*n){
        risk.over.boin[count] <- 1
      }
    }else if(true.mtd==K){
      if(sum(dose.treated[1:(K-1)])*cohortsize>risk.cutoff*n){
        risk.under.boin[count] <- 1
      }
    }else{
      if(sum(dose.treated[1:(true.mtd-1)])*cohortsize>risk.cutoff*n){
        risk.under.boin[count] <- 1
      }else if(sum(dose.treated[(true.mtd+1):K])*cohortsize>risk.cutoff*n){
        risk.over.boin[count] <- 1
      }
    }
  }
  
  ## Prepare output
  temp <- data.frame(matrix(nrow = 8,ncol = K+6))
  temp[1,] <- c(paste("Scenario",s),"Pr(toxicity)",skeleton,NA,NA,NA,NA)
  temp[2,] <- c("PBF","Selection(%)",round(c(cont(mtd.pbf,n.level = K)/n.sim,
                        mean(risk.over.pbf),mean(risk.under.pbf)),digits = 3)*100,mean(dose.es),mean(dose.de))
  temp[3,] <- c(NA,"# Patients",round(num.p.pbf*cohortsize/n.sim,digits = 1),NA,NA,NA,NA)
  temp[4,] <- c("BOIN","Selection(%)",round(c(cont(mtd.boin,n.level = K)/n.sim,
                        mean(risk.over.boin),mean(risk.under.boin)),digits = 3)*100,mean(dose.es.boin),mean(dose.de.boin))
  temp[5,] <- c(NA,"# Patients",round(num.p.boin*cohortsize/n.sim,digits = 1),NA,NA,NA,NA)
  temp[6,] <- c("PBF-next","Selection(%)",round(cont(mtd.next.pbf,n.level = K)/n.sim,digits = 3)*100,NA,NA,NA,NA)
  temp[7,] <- c("BOIN-next","Selection(%)",round(cont(mtd.next.boin,n.level = K)/n.sim,digits = 3)*100,NA,NA,NA,NA)
  tab <- rbind(tab,temp)
}

write_xlsx(tab,"D:/study/Dissertation/PrBF design/Figure/Numbers1.xlsx",
                      col_names = FALSE)


