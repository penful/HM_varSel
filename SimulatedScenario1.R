#########################################################
#### Code provided with the scientific paper titled: 
### "Variable selection for hidden Markov models with continuous variables and missing data" 
### Authors: F. Pennoni, F. Bartolucci, F. Pandolfi 
### Journal of Classification (2023)
#########################################################

# Required packages
rm(list=ls())
require(clustvarsel)
require(LMest)
library(mvtnorm)
require(Formula)
require(MASS)
require(mix)

# load the required fuctions written for the proposal
source("regress_miss.R")         # Function to fit the multivariate regression model
source("lmestContMISS.R")        # Function to estimate the HM model with missing data
source("functions.R")            # Internal function
source("lmbasic.cont.MISS.R")    # Internal function
source("complk_cont_miss.R")     # Internal function
source("drawHMBasicCont.R")      # Function to simulate data
source("compute_BIC.R")          # Function to compute BIC
source("item_selection.R")       # Function to select responses
source("count_eq.R")             # Internal function
source("forward_regress_miss.R") # Function to apply forward model selection

#### simulate data ####

# B to be changed for number of simulated samples B = 100
B <- 1
# n to be changed number of individuals
n <- 25   # number of individuals to increase  (250,500)
TT <- 5 # time occasions (5,10)
k <- 2  # number of hidden states (2,3)
J <- 5 # total number of variables 
r <- 2 # clustering variables (2,4)
pmiss <- 0.05 # missing proportion (0.1,0.25)
nmiss <- pmiss*(n*TT)
piv <- rep(1/k,k) # initial probabilities
modBasic <- 1 # time homogenous transitions
Kmax <- 5
ind <- TRUE
if(k==2){
  Pi <- matrix(c(0.8,0.2,0.2,0.8), k, k,byrow=TRUE)
  Pi <- array(Pi, c(k, k, TT))
  Pi[,,1] <- 0
  if(r==2) Mu <- matrix(c(0,0,4,0),r,k)
  else if(r==4) Mu <- matrix(c(0,0,0,0,4,0,4,0), r, k)
}else if(k==3){
  Pi <- matrix(c(0.80,0.10,0.10,0.10,0.80,0.10,0.10,0.10,0.80), k, k,byrow=TRUE)
  Pi <- array(Pi, c(k, k, TT))
  Pi[,,1] <- 0
  if(r==2) Mu <- matrix(c(0,0,4,0,4,2),r,k)
  else if(r==4)  Mu <- matrix(c(0,0,0,0,4,0,4,0,4,2,4,2), r, k)
}
Si <- matrix(0.5,r,r)
diag(Si) <- 1

# file generated with the results along with another partial file of the results # 
filename <- sprintf("Results_n%g_k%g_TT%g_r%g_J%g_pmiss%g_ind%s.RData",n,k,TT,r,J,pmiss,ind)
sim = res <- vector("list",B)
ARI = timeE<- rep(0,B)

for(b in 1:B){
  try({
    print(paste("sample:", b))
   set.seed(b+15200)
    sim[[b]] <- drawHMBasicCont(piv, Pi, Mu, Si, n, 
                           format = "long")
    
    Yc <- sim[[b]]$Y
    
    if(ind){
      a <- seq(-2,2,length.out=J-r)
      Yind <- rmvnorm(n*TT,a,diag(J-r))
      Y.final <- cbind(Yc,Yind)
    } 
    else{
      a <- seq(-2,2,length.out=J-r)
      be <- matrix(c(-1,1),r,J-r)
      Yc2 <- as.matrix(Yc[,3:(r+2)])
      Yreg <- rep(1,n*TT)%x%t(a) +Yc2%*%be+rmvnorm(n*TT,rep(0,J-r),diag(J-r))
      Y.final <- cbind(Yc,Yreg)
    }
  
#### Simulate missing values under the Missing-at-Random assumption ####
    if(pmiss>0){
      Y.final2 <- Y.final
      for(j in 3:(J+2)) Y.final2[runif(n*TT)<pmiss,j] = NA
    }else{
      Y.final2 <- Y.final
    }


#### Perform model and variable selection ####
    Y <- data.frame(Y.final2)
    time <- proc.time()
    res[[b]] <- item_selection(Y,modBasic=modBasic)
    time <- proc.time()-time
    timeE[b] <- time[3]	
    
    class <- apply(res[[b]]$mod_HM_now$V,c(1,3),which.max)
    ARI[b] <- adjustedRandIndex(sim[[b]]$U,class)
  }) 
  if(b/10==floor(b/10)) save.image(filename)
}

save.image(filename)

