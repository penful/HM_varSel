item_selection <- function(Y,modBasic=1,Kmax=5,tol=10^-10){

#---- preliminaries ----
  time = proc.time()
  id = colnames(Y)[1]
  tt = colnames(Y)[2]
  J = ncol(Y)-2
  n = length(unique(Y[,1]))

#---- find starting step ----
  cat("*** FIND INITIAL ITEM TO INCLUDE ***\n")
  bicv = rep(0,J); bic_best = Inf; out_best=NULL
  Yimp = NULL
  for(j in 1:J){

    out = compute_BIC(Y,j,(1:J)[-j],1:Kmax,Kmax=Kmax,Yimp=Yimp,tol=tol,modBasic = modBasic) 
    Yimp = out$Yimp
    bicv[j] = out$bic_tot
    if(out$bic_tot<bic_best){
      bic_best = out$bic_tot
      out_best = out
      mod_HM_best = out$out_HM
      k_best = out$k
    }
    cat("*** it_glob = ",0,", it = ",j,", j = ",j,", bic = ",out$bic_tot,
        ", bic_best = ",bic_best,"\n",sep="")
  }
  items_now = which.min(bicv)
  mod_HM_now = mod_HM_best
  bic_now = bic_best
  k_now = k_best

#---- iterate until convergence ----
  cat("*** STEPWISE SELECTION OF ITEMS ***\n")
  flag = TRUE
  it_glob = 0
  while(flag){
    cat("items_now =",items_now,", k_now =",k_now,"\n")
# inclusion step
    it_glob = it_glob+1	
    print("inclusion step")
    items_no = setdiff(1:J,items_now) 
    bic_diffv = rep(0,length(items_no))
    it = 0
    bic_diff_best = Inf
    for(j in items_no){
      it = it+1
      items_tmp = sort(c(items_now,j))
      if(k_now==1){
        kv = 1:2
      }else if(k_now==Kmax){
        kv = c(Kmax-1,Kmax)
      }else{
        kv = c(k_now-1,k_now,k_now+1)
      }
      
      out = compute_BIC(Y,items_tmp,(1:J)[-items_tmp],kv,Yimp=Yimp,tol=tol,modBasic = modBasic)

      Yimp = out$Yimp
      Y_cov <- Yimp[,items_now,drop=FALSE]
      Y_cov1 = Y_cov[,1]
      if(ncol(Y_cov)>1) for(j1 in 2:ncol(Y_cov)){
        Tmp = cbind(1,Y_cov1,Y_cov[,j1])
        if(rcond(t(Tmp)%*%Tmp)>10^-15) Y_cov1 = cbind(Y_cov1,Y_cov[,j1])
      }
      mod_reg <- forward_regress_miss(as.matrix(Y[,items_no+2]),as.matrix(Y_cov1),it,n1=NULL,tol=tol) 
      bic_diffv[it] = out$bic_tot-(mod_HM_now$bic+mod_reg$bic)
      if(bic_diffv[it]<bic_diff_best){
        bic_diff_best = bic_diffv[it]
        mod_HM_best = out$out_HM
        k_best = out$k
        out_best = out
      }
      cat("it_glob = ",it_glob,", it = ",it,", items_now = ",paste(items_now,collapse=",",sep=""),
          ", j = ",j,", k_now = ",k_now,", bic_now = ",bic_now,", bic_best = ",bic_best,"\n",sep="")
    }

    if(bic_diff_best<0){
      item_add = items_no[which.min(bic_diffv)]
      items_now = c(items_now,item_add)
      mod_HM_now = mod_HM_best
      k_now = k_best
      bic_now = out_best$bic_tot
      flag = TRUE
      print("==> include item")
      print(item_add)
      
	  save.image("preliminary_new.RData")	
    }else flag = FALSE
    cat("items_now =",items_now,", k_now =",k_now,"\n")
    
# exclusion step
    print("exclusion step")
    
    bicv = rep(0,length(items_now))
    it = 0
    bic_best = Inf
    for(j in items_now){
      it = it+1
      items_tmp = sort(setdiff(items_now,j))
      if(k_now==1) kv = 1:2 else kv = c(k_now-1,k_now,k_now+1)
      
      out = compute_BIC(Y,items_tmp,(1:J)[-items_tmp],kv,Yimp=Yimp,tol=tol,modBasic = modBasic)
      Yimp = out$Yimp
      
      Y_cov <- Yimp[,items_tmp,drop=FALSE]
      Y_cov1 = Y_cov[,1]
      if(ncol(Y_cov)>1) for(j1 in 2:ncol(Y_cov)){
        Tmp = cbind(1,Y_cov1,Y_cov[,j1])
        if(rcond(t(Tmp)%*%Tmp)>10^-15) Y_cov1 = cbind(Y_cov1,Y_cov[,j1])
      }
      mod_reg <- forward_regress_miss(as.matrix(Y[,(1:J)[-items_tmp]+2]),as.matrix(Y_cov1),
                                      which((1:J)[-items_tmp]==j),n1=NULL,tol=tol)  
      bicv[it] <- out$bic_HM+mod_reg$bic
      bic_tmp <- out$bic_HM+mod_reg$bic
      if(bic_tmp<bic_best){
        bic_best = bic_tmp
        mod_HM_best = out$out_HM
        k_best = out$k
      }
      cat("it_glob = ",it_glob,", it = ",it,", items_now = ",paste(items_now,collapse=",",sep=""),
          ", j = ",-j,", k_now = ",k_now,", bic_now = ",bic_now,", bic_best = ",bic_best,"\n",sep="")
    }

    if(bic_best<bic_now){
      item_drop = items_now[which.min(bicv)]
      items_now = setdiff(items_now,item_drop)
      mod_HM_now = mod_HM_best
      k_now = k_best
      bic_now = bic_best
      flag = TRUE
      print("==> exclude item")
      print(item_drop)
 #  	  save.image("preliminary_new.RData")	
    }

#save.image("item_sel_prel.RData")
  }
  time = proc.time()-time
  out = list(mod_HM_now=mod_HM_now,items=items_now,bic=bic_now,k=k_now,time=time,Yimp=Yimp)

}