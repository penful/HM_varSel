forward_regress_miss <- function(Y,X,ind,n1=NULL,tol=10^-10){

# stepwise regression for Y with the covariates X
# ind is the index of the response variable on which it makes model selection

#---- preliminaries ----
  n = nrow(Y)
  ncx = ncol(X)
  ncy = ncol(Y)

#---- initial model ----
  XX = vector("list",ncy)
  for(j in (1:ncy)[-ind]) XX[[j]] = X
  XX[[ind]] = matrix(0,n,0)
  mod = regress_miss(Y,XX,n1=n1,disp = TRUE,reltol=tol)
  covin = NULL
  print(mod$bic)
  print("-----")

# ---- try tro introduce a new covariate ----
  flag = TRUE
  while(flag){
    covin0 = covin
    flag = FALSE
    tmp = 1:ncx
    if(!is.null(covin0)) tmp = tmp[-covin]
    for(h in tmp){
      XX[[ind]] = as.matrix(X[,c(covin0,h)])
      mod1 = regress_miss(Y,XX,n1=n1,disp = TRUE,reltol=tol)
      print(c(covin0,h,mod1$bic))
      if(mod1$bic<mod$bic){
        flag = TRUE
        mod = mod1
        covin = c(covin0,h)
      }
    }
    if(flag) cat("* add",covin[length(covin)],"*\n")
    print("-----")
  }

#---- output ----
  out = mod
  out$covin = covin
  return(out)

}