drawHMBasicCont <-
function(piv,Pi,Mu,Si,n,format = c("long","matrices")){

#        [Y,yv] = draw_lm_basic(piv,Pi,Mu,Si,n)
#
# Draw a sample of size n from a Basic Latent Markov model for continuous data with parameter piv, Pi, Mu and Si

  
  # Preliminaries
  format <- match.arg(format, choices = eval(formals(drawLMbasiccont)$format))
  

# Preliminaries
	if(is.vector(Mu)){
    	r = 1
    	k = length(Mu)
    	Mu = matrix(Mu,r,k)
    }else{
    	r = nrow(Mu)
    	k = ncol(Mu)
    }
	TT = dim(Pi)[3]
	if(r==1) Si = matrix(Si,r,r)
# For each subject
    Y = array(0,c(n,TT,r))
    U = matrix(0,n,TT)
    cat("------------|\n")
    cat(" sample unit|\n")
    cat("------------|\n")
    for(i in 1:n){
      if(i/1000==floor(i/1000)) cat(sprintf("%11g",i),"\n",sep=" | ")
      u = k+1-sum(runif(1)<cumsum(piv))
      Y[i,1,] = rmvnorm(1,Mu[,u],Si)
      U[i,1] = u
      for(t in 2:TT){
        u = k+1-sum(runif(1)<cumsum(Pi[u,,t]))
        Y[i,t,] = rmvnorm(1,Mu[,u],Si)
        U[i,t] = u
      }
    }
    if(format == "long")       Y <- matrices2long(Y = Y)
    
    cat("------------|\n")
    	out = list(Y=Y,U=U,piv = piv, Pi = Pi, Mu = Mu, Si = Si, n = n)
}
