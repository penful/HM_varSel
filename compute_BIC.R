compute_BIC <- function(Y,items_HM,items_reg,k,modBasic=1,Kmax=5,Yimp=NULL,tol=10^-10){

#---- preliminaries ----
  n = length(unique(Y[,1]))

#---- imputation if necessary ----
  miss <- any(is.na(Y))
  if(miss){
    if(is.null(Yimp)){
      cat("*** DATA IMPUTATION ***\n")
      formula <- lmestFormula(data=data.frame(Y[,-(1:2)]),response=1:(ncol(Y)-2))
      mod0 <- lmestContMISS(index=names(Y)[1:2], k = 1, data = Y,
                            responsesFormula = formula$responsesFormula,modBasic=modBasic,tol=tol)
      k0 = 1
      Yimp <- mod0$Yimp
      Ydim <- dim(Yimp)
      Yimp <- matrix(aperm(Yimp,c(2,1,3)),Ydim[1]*Ydim[2],Ydim[3])
      flag = TRUE
      cat("*** k =",1,"BIC = ",mod0$bic,"\n")
      while(flag){
        flag = FALSE
        out0 = kmeans(Yimp,k0+1,nstart = 200)
        U = t(matrix(out0$cluster,Ydim[2],Ydim[1]))
        Mu = t(out0$centers)
        Si = matrix(0,Ydim[3],Ydim[3])
        for(u in 1:(k0+1)){
          Tmp = Yimp[out0$cluster==u,]-rep(1,sum(out0$cluster==u))%o%Mu[,u]
          Si = Si+t(Tmp)%*%Tmp/(Ydim[1]*Ydim[2])
        }
        piv = rep(0,(k0+1))
        for(u in 1:(k0+1)) piv[u] = sum(U[,1]==u)/Ydim[1]
        Pi = matrix(0,k0+1,k0+1)
        for(u in 1:(k0+1)) for(v in 1:(k0+1)) Pi[u,v] = sum(U[,1:(Ydim[2]-1)]==u & U[,2:Ydim[2]]==v)
        Pi = (1/rowSums(Pi))*Pi
        Pi = array(Pi,c(k0+1,k0+1,Ydim[2])); Pi[,,1] = 0
        mod1 <- lmestContMISS(index=names(Y)[1:2], k = k0+1, data = Y,
                              responsesFormula = formula$responsesFormula,
                              modBasic=modBasic,tol=tol,
                              parInit = list(piv=piv,Pi=Pi,Mu=Mu,Si=Si),start = 2)
        cat("*** k =",k0+1,"BIC = ",mod1$bic,"\n")
        if(mod1$bic<mod0$bic){
          mod0 = mod1
          k0 = k0+1
          if(k0<Kmax) flag = TRUE
          Yimp <- mod0$Yimp
          Ydim <- dim(Yimp)
          Yimp <- matrix(aperm(Yimp,c(2,1,3)),Ydim[1]*Ydim[2],Ydim[3])
        }
      }
    }
  }else{
    Yimp <- Y[,-(1:2)]
  }

# fit HM model
  Y_HM <- Y[,c(1,2,items_HM+2)]
  formula <- lmestFormula(data=data.frame(Y[,-c(1,2)]),response=items_HM)
  k0 = k[1]
  nt = nrow(Yimp)
  TT = nt/n
  r_HM = length(items_HM)
  if(k[1]==1){
    out_HM <- lmestContMISS(index=names(Y)[1:2], k = k0, data = Y_HM,
                            responsesFormula = formula$responsesFormula,modBasic=modBasic,tol=tol)
  }else{
    Yimp1 = as.matrix(Yimp[,items_HM])
    out1 = kmeans(Yimp1,k0,nstart = 200)
    U = t(matrix(out1$cluster,TT,n))
    Mu = t(out1$centers)
    Si = matrix(0,r_HM,r_HM)
    for(u in 1:k0){
      Tmp = Yimp1[out1$cluster==u,]-rep(1,sum(out1$cluster==u))%o%Mu[,u]
      Si = Si+t(Tmp)%*%Tmp/nt
    }
    piv = rep(0,k0)
    for(u in 1:k0) piv[u] = sum(U[,1]==u)/n
    Pi = matrix(0,k0,k0)
    for(u in 1:k0) for(v in 1:k0) Pi[u,v] = sum(U[,1:(TT-1)]==u & U[,2:TT]==v)
    Pi = (1/rowSums(Pi))*Pi
    Pi = array(Pi,c(k0,k0,TT)); Pi[,,1] = 0
    out_HM <- lmestContMISS(index=names(Y)[1:2], k = k0, data = Y_HM,
                             responsesFormula = formula$responsesFormula,modBasic=modBasic,tol=tol,
                             parInit = list(piv=piv,Pi=Pi,Mu=Mu,Si=Si),start = 2)
  }
  if(length(k)>1) for(k0 in k[-1]){
    if(k0==1){
      out_HM1 <- lmestContMISS(index=names(Y)[1:2], k = k0, data = Y_HM,
                               responsesFormula = formula$responsesFormula,modBasic=modBasic,tol=tol)
    }else{
      Yimp1 = as.matrix(Yimp[,items_HM])
      out1 = kmeans(Yimp1,k0,nstart = 200)
      U = t(matrix(out1$cluster,TT,n))
      Mu = t(out1$centers)
      Si = matrix(0,r_HM,r_HM)
      for(u in 1:k0){
        Tmp = Yimp1[out1$cluster==u,]-rep(1,sum(out1$cluster==u))%o%Mu[,u]
        Si = Si+t(Tmp)%*%Tmp/nt
      }
      piv = rep(0,k0)
      for(u in 1:k0) piv[u] = sum(U[,1]==u)/n
      Pi = matrix(0,k0,k0)
      for(u in 1:k0) for(v in 1:k0) Pi[u,v] = sum(U[,1:(TT-1)]==u & U[,2:TT]==v)
      Pi = (1/rowSums(Pi))*Pi
      Pi = array(Pi,c(k0,k0,TT)); Pi[,,1] = 0
      out_HM1 <- lmestContMISS(index=names(Y)[1:2], k = k0, data = Y_HM,
                               responsesFormula = formula$responsesFormula,modBasic=modBasic,tol=tol,
                               parInit = list(piv=piv,Pi=Pi,Mu=Mu,Si=Si),start = 2)
    }

    if(out_HM1$bic<out_HM$bic) out_HM = out_HM1
  }

# fit multiple linear regression
  miss <- any(is.na(Y))
  Y_reg <- Y[,items_reg+2]
  Y_cov <- Yimp[,items_HM,drop=FALSE]
  Y_cov1 = Y_cov[,1]
  if(ncol(Y_cov)>1) for(j in 2:ncol(Y_cov)){
    Tmp = cbind(1,Y_cov1,Y_cov[,j])
    if(rcond(t(Tmp)%*%Tmp)>10^-15) Y_cov1 = cbind(Y_cov1,Y_cov[,j])
  }
  out_reg <- regress_miss(as.matrix(Y_reg),as.matrix(Y_cov1),n1=NULL,reltol=tol,disp=TRUE)
# results
  lk_HM = out_HM$lk; lk_reg = out_reg$lk
  np_HM = out_HM$np; np_reg = out_reg$np
  bic_HM = out_HM$bic; bic_reg = out_reg$bic
  if(length(lk_reg)==0) lk_tot = lk_HM
  else lk_tot = lk_HM+lk_reg
  np_tot = np_HM+np_reg
  if(length(bic_reg)==0) bic_tot = bic_HM
  else bic_tot = bic_HM+bic_reg

  res = data.frame(HM=c(lk_HM,np_HM,bic_HM),reg=c(lk_reg,np_reg,bic_reg),
                   tot=c(lk_tot,np_tot,bic_tot))
  rownames(res) = c("lk","np","BIC")
  k = out_HM$k
  out = list(lk_HM=lk_HM,lk_reg=lk_reg,np_HM=np_HM,np_reg=np_reg,
             bic_HM=bic_HM,bic_reg=bic_reg,lk_tot=lk_tot,np_tot=np_tot,
             out_HM=out_HM,out_reg=out_reg,
             bic_tot=bic_tot,res=res,k=k,Y_reg=Y_reg,Y_cov=Y_cov,Yimp=Yimp)
  return(out)

}
