regress_miss <- function(Y,X,disp=FALSE,itmax=1000,reltol=10^-10,n1=NULL,Yimp1=NULL){

# fit regression model (also multivariate) under MAR on y
# intercept is automatically included
#
# INPUT:
# Y = matrix of responses
# X = design matrix or list of design matrices of each response

# Preliminaries
  ldnorm1 <-function(y,mu,si2) lf = -(y-mu)^2/si2/2-log(2*pi*si2)/2
  ldmvnorm1 <-function(y,mu,Si) lf = -c((y-mu)%*%solve(Si)%*%(y-mu))/2-log(det(2*pi*Si))/2
  Y = as.matrix(Y)
  n = nrow(Y)
  if(is.list(X)){
    for(j in 1:ncol(Y)) X[[j]] = cbind(1,X[[j]])
  }else{
    X = cbind(1,X)
  }
  if(is.null(n1)) n1 = n

### multiple regression
  if(ncol(Y)==1){
    if(is.list(X)) X = X[[1]]
# fit model for non-missing data
    ind = !is.na(Y)
    modlm = lm(Y~-1+X)
    be = modlm$coefficients
    lk = logLik(modlm)[1]
# compute likelihood
    Yimp = Y
    if(any(!ind)) Yimp[!ind] = X[!ind,,drop=FALSE]%*%be
    np = length(be)+1
    bic = -2*lk + log(n1)*np
    si2 = summary(modlm)$sigma^2
# output
    out = list(be=be,si2=si2,lk=lk,np=np,bic=bic,yimp=Yimp)
  }else{
### multivariate regression
    miss = any(is.na(Y))
    ncy = ncol(Y)
    if(is.list(X)){
      ncx = rep(0,ncy)
      for(j in 1:ncy) ncx[j] = ncol(X[[j]])
    }else{
      ncx = ncol(X)
    }
    if(is.list(X)){
      XX = array(0,c(ncy,sum(ncx),n))
      ind = 0
      for(j in 1:ncy){
        ind = max(ind)+(1:ncx[j])
        XX[j,ind, ] = t(X[[j]])  
      }
    }

## with missing data
    if(miss){
      R = (!is.na(Y))
# starting values
      if(is.null(Yimp1)) Yimp2 = Y else Yimp2 = Yimp1
      B = matrix(0,ncx,ncy)
      if(is.list(X)){
        B = vector("list",ncy)
        Mu = matrix(0,n,ncy)
        for(j in 1:ncy){
          if(nrow(X[[j]])==1) out = regress_miss(Yimp2[,j],matrix(0,n,0)) else out = regress_miss(Yimp2[,j],X[[j]][,-1])
          B[[j]] = out$be
          Mu[,j] = X[[j]]%*%B[[j]]
        }
      }else{
        for(j in 1:ncy){
          out = regress_miss(Yimp2[,j],X[,-1])
          B[,j] = out$be
        }
        Mu = X%*%B
      }
      if(is.null(Yimp1)){
        Yimp = Y
        for(i in 1:n){
          indo = R[i,]
          if(sum(indo)<ncy) Yimp[i,!indo] = Mu[i,!indo]
        }
      }else{
        Yimp = Yimp1
      }
      E = Yimp-Mu
      Si = (t(E)%*%E)/n
# compute log-likelihood
      lk = 0
      for(i in 1:n){
        indo = R[i,]
        if(sum(indo)==1){
          lk = lk+ldnorm1(Y[i,indo],Mu[i,indo],Si[indo,indo])
        }else if(sum(indo)>1){
          lk = lk+ldmvnorm1(Y[i,indo],Mu[i,indo],Si[indo,indo])
        }
      }
      it = 0; lko = lk
      if(disp){
        cat("------------|-------------|-------------|\n")
        cat("    step    |     lk      |    lk-lko   |\n")
        cat("------------|-------------|-------------|\n")
        cat(sprintf("%11g",c(0,lk)),"\n",sep = " | ")
      }
      while(((lk-lko)/abs(lko)>reltol | it==0) & it<itmax){
        # t0 = proc.time()
        it = it+1
# E-step
        Yimp = Y; Vc = matrix(0,ncy,ncy)
        sing = 0
        for(i in 1:n){
          indo = R[i,]
          if(sum(indo)==0){
            Yimp[i,] = Mu[i,]
            Vc = Vc+Si
          }else if(sum(indo)<ncy){
            iSi = try(solve(Si[indo,indo]),silent=TRUE)
            if(inherits(iSi,"try-error")){
              sing = sing+1
              iSi = ginv(Si[indo,indo])
            }
            Tmp = Si[!indo,indo]%*%iSi
            Yimp[i,!indo] = Mu[i,!indo]+Tmp%*%(Y[i,indo]-Mu[i,indo])
            Vc[!indo,!indo] = Vc[!indo,!indo]+Si[!indo,!indo]-Tmp%*%Si[indo,!indo]
          }
        }
        if(sing>0) cat("***",sing,"singular matrices ***\n")
        # print(c(1,proc.time()-t0))
# M-step
       if(is.list(X)){
         iSi = try(solve(Si),silent=TRUE)
         if(inherits(iSi,"try-error")) iSi = ginv(Si)
         
         NUM = DEN = 0
         for(i in 1:n){
           Tmp = t(XX[,,i])%*%iSi
           NUM = NUM + Tmp%*%Yimp[i,]
           DEN = DEN + Tmp%*%XX[,,i]
         }
         iDEN = try(solve(DEN),silent=TRUE)
         if(inherits(iDEN,"try-error")) iDEN = ginv(DEN)
         be = iDEN%*%NUM
         ind = 0
         for(j in 1:ncy){
           ind = max(ind)+(1:ncx[j])
           B[[j]] = be[ind]
           Mu[,j] = X[[j]]%*%B[[j]]
         }
       }else{
         M = t(X)%*%X
         B = solve(M,t(X))%*%Yimp
         Mu = X%*%B
       }
       E = Yimp-Mu
       Si = (Vc+t(E)%*%E)/n
       # print(c(2,proc.time()-t0))
# compute log-likelihood
        lko = lk; lk = 0
        for(i in 1:n){
          indo = R[i,]
          if(sum(indo)==1){
            lk = lk+ldnorm1(Y[i,indo],Mu[i,indo],Si[indo,indo])
          }else if(sum(indo)>1){
            lk = lk+ldmvnorm1(Y[i,indo],Mu[i,indo],Si[indo,indo])
          }
        }
        # print(c(it,lk,lk-lko))
        # print(c(3,proc.time()-t0))
         if(disp & it%%100==0) cat(sprintf("%11g",c(it,lk,lk-lko)),"\n",sep = " | ")
      }
      if(disp & it%%100>0) cat(sprintf("%11g",c(it,lk,lk-lko)),"\n",sep = " | ")
      if(disp) cat("------------|-------------|-------------|\n")
    }else{
## without missing data
      if(is.list(X)){
        B = vector("list",ncy)
        Mu = matrix(0,n,ncy)
        for(j in 1:ncy){
          tmp = try(solve(t(X[[j]])%*%X[[j]],t(X[[j]])),silent=TRUE)
          if(inherits(tmp,"try-error")) tmp = ginv(t(X[[j]])%*%X[[j]])%*%t(X[[j]])
          B[[j]] = c(tmp%*%Y[,j])
          Mu[,j] = X[[j]]%*%B[[j]]
        }
        E = Y-Mu
        Si = (t(E)%*%E)/n
        lk = 0
        for(i in 1:n) lk = lk+dmvnorm(Y[i,],Mu[i,],Si,log=TRUE)
        lko = lk; it = 0
        while(abs(lk-lko)/abs(lko)>10^-8|it == 0){
          lko = lk; it = it+1
          iSi = try(solve(Si),silent=TRUE)
          if(inherits(iSi,"try-error")) iSi = ginv(Si)
          NUM = DEN = 0
          for(i in 1:n){
            Tmp = t(XX[,,i])%*%iSi
            NUM = NUM + Tmp%*%Y[i,]
            DEN = DEN + Tmp%*%XX[,,i]
          }
          iDEN = try(solve(DEN),silent=TRUE)
          if(inherits(iDEN,"try-error")) iDEN = ginv(DEN)
          be = iDEN%*%NUM
         
          ind = 0
          for(j in 1:ncy){
            ind = max(ind)+(1:ncx[j])
            B[[j]] = be[ind]
            Mu[,j] = X[[j]]%*%B[[j]]
          }
          E = Y-Mu
          Si = (t(E)%*%E)/n
          lk = 0
          for(i in 1:n) lk = lk+dmvnorm(Y[i,],Mu[i,],Si,log=TRUE)
        }
      }else{
        B = solve(t(X)%*%X,t(X))%*%Y
        Mu = X%*%B
        E = Y-Mu
        Si = (t(E)%*%E)/n
        lk = 0
        for(i in 1:n) lk = lk+dmvnorm(Y[i,],Mu[i,],Si,log=TRUE)
      }
      Yimp = Y
    }
    if(is.list(X)){
      np = sum(ncx)+ncy*(ncy+1)/2
    }else{
      np = ncy*ncx+ncy*(ncy+1)/2
    }
    bic = -2*lk+log(n1)*np
    out = list(B=B,Si=Si,lk=lk,np=np,bic=bic,Yimp=Yimp)
  }

# output
  return(out)

}