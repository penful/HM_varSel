rm(list=ls())
load("Results_n25_k2_TT5_r2_J5_pmiss0.05_indTRUE.RData")
items.sel <- matrix(0,B,J)
k.sel <- rep(0,B)
sel = selin <- rep(0,B)
for(b in 1:B){
  items.sel[b,1:length(res[[b]]$items)] <- res[[b]]$items
  k.sel[b] <- res[[b]]$k
 if(setequal(items.sel[b,1:length(res[[b]]$items)],1:r)) sel[b] = 1
 if(sum((1:r)%in%items.sel[b,])==r) selin[b] = 1
}
print("%%%%%%%%%%%%%%%%%%%%%%%%")
print("Results starting from a single item")
print("Frequency of true k")
print(table(k.sel)/B)
print("Frequency of true items")
print(sum(sel)/B)
print("Frequency of samples containing the items")
print(sum(selin)/B)
print("Average computing time")
print(mean(timeE))
print("Average ARI")
print(mean(ARI))