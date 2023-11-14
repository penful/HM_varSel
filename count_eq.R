count_eq <- function(A,b){
  
  eq = 0; cont = 0
  if(is.vector(A)) A = t(A)
  n = nrow(A)
  for(i in 1:n){
    a = A[i,]
    a = a[a>0]; b = b[b>0]
    eq = eq + (all(a%in%b) & all(b%in%a))/n
    cont = cont + all(b%in%a)/n
  }
  res = c(eq=eq,cont=cont)
  return(res)

}