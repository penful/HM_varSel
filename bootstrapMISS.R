bootstrapMISS <- function(Y,Mu,Si,piv,Pi,modBasic,B=100){

# preliminaries
	mMu = mSi = 0
	m2Mu = m2Si  = 0
	mpiv = mPi = 0
	m2piv = m2Pi = 0
	
  if(is.vector(Mu)){
    r =1
    k = length(Mu)
  }else{
    r = nrow(Mu)
    k = ncol(Mu)
  }
	units <- unique(Y$id)
	n <- length(units)

  for (b in 1:B) {
    cat("non-parametric boostrap sample n. ",b,"\n")
   ind = sample(n,n,replace=T)
    unitsB <- units[ind]
    Yb = Y[Y$id==unitsB[1],]
    for(i in 2:n){
      Yb = rbind(Yb,Y[Y$id==unitsB[i],])
    }
      
    formula <- lmestFormula(data=Yb[,-c(1,2)],response=1:ncol(Yb[,-c(1,2)]))
    out <- lmestContMISS(index=c("id","time"),
                          k = k, data = Yb,
                          responsesFormula = formula$responsesFormula, modBasic=modBasic,
                         start=2,parInit = list(piv=piv,Pi=Pi,Mu=Mu,Si=Si))
    mMu = mMu + out$Mu/B
    mSi = mSi + out$Si/B
    mpiv = mpiv+out$piv/B
    mPi = mPi+out$Pi/B
    
    m2Mu = m2Mu + out$Mu^2/B
    m2Si = m2Si + out$Si^2/B
    m2piv = m2piv+out$piv^2/B
    m2Pi = m2Pi+out$Pi^2/B
    if(b/10==floor(b/10)) save.image("boot_tmp.RData")
  }
  seMu = sqrt(m2Mu - mMu^2)
  seSi = sqrt(m2Si - mSi^2)
  sepiv = sqrt(m2piv-mpiv^2) 
  sePi = sqrt(m2Pi-mPi^2)
  
  out = list(mMu = mMu, mSi = mSi, mpiv = mpiv, mPi = mPi,
                 seMu = seMu, seSi = seSi, sepiv = sepiv, sePi = sePi)
  
}