# Read results of the estimated model 

rm(list=ls())

load("example_trial.Rdata")

#######  Descriptive  Plots of the variables ################################
####  Plots and summary of the data ####
###############################################

require(ggplot2)
ggplot(dt1, 
       aes(x = time, 
           y = as.factor(id), 
           fill =  Y1))+ 
  geom_tile(colour = 'gray98') +
  scale_fill_gradient(na.value = 'black',
                      #na.value = 'white', 
                      low = "yellow", high = "darkgreen") +
  scale_x_discrete(expand = c(0, 0)) +    
  scale_y_discrete(expand = c(0, 0), #labels = NULL, 
                   breaks = NULL) +
  labs(x = " ", y = "region") + 
  theme(
    plot.background=element_rect(fill="black", color=NA),
    panel.background=element_rect(fill="white", color=NA),
    #axis.text.y = element_text(data1$Country_Code),
    # axis.text.x = element_blank(),
    #  legend.title=element_blank(),
    legend.key = element_rect(colour = "black", fill = 'NA')
  )


#######  Model results ################################
####  Estimated Means  ####
###############################################
mod.final <- out1$mod_HM_now

# fix the estimated parameters

k.final <- mod.final$k
Mu.final <- mod.final$Mu
piv.final <- mod.final$piv
Pi.final <- mod.final$Pi
Si.final <- mod.final$Si
Y.final <- mod.final$data
dim(Y.final)

# selected items
names(Y.final)

#### Description of selected indicators ####

MIN = apply(Y.final[,3:7],2,min,na.rm=TRUE); MIN
Q1 = apply(Y.final[,3:7],2,quantile,0.25,na.rm=TRUE)
ME = apply(Y.final[,3:7],2,median,na.rm=TRUE)
Q3 = apply(Y.final[,3:7],2,quantile,0.75,na.rm=TRUE)
MX = apply(Y.final[,3:7],2,max,na.rm=TRUE)
SD = apply(Y.final[,3:7],2,sd,na.rm=TRUE)
A <- cbind(MIN,Q1, ME, Q3, MX, SD)
round(A,2)

# Estimated cluster means
round(mod.final$Mu,2)
# Order for  increasing means of item 2
indC <- order(mod.final$Mu[2,])
indC
Munew = mod.final$Mu[,indC]
round(Munew,3)


# Heat map once the values are standardized by row
norm01 <- function(x) (x - mean(x))/sd(x)
M <- apply(t(Munew), 2, norm01)
require("RColorBrewer")
heatmap(t(M), Rowv = NA, Colv = NA, 
        col= colorRampPalette(brewer.pal(8, "BuGn"))(25),
        scale = "none", 
        labCol = paste("State", 1:k.final), cexCol = 1.1)


######  Table ################################
####  Estimated Variances  ####
###############################################

SS <- round(mod.final$Si,3)
require(ggm)
SSa<-round(parcor(SS),3)
require(psych)
AA<-lowerUpper(SS, SSa)
# variances on the diagonal
diag(AA)<-round(diag(SS),3)
round(AA, 3)
AA

#######  Table ################################
####  Initial and transition probabilities  ####
###############################################
# Initial probabilities 
print(round(mod.final$piv[indC],3))
# Transition probabilities 
round(mod.final$Pi[, ,5],3)

#######  Bootstrap s.e. ################################
####  Report bootstrap s.e.  ####
###############################################

# Perform non-parametric bootstrap
source("bootstrapMISS.R")


# fix the estimated parameters

head(Y.final)
# The first two columns must be "id" and "time" 

# fix the number of bootstrap replicates: note that B should be 200 or more  
B <- 5
resboot <- bootstrapMISS(Y.final,Mu.final,Si.final,
                         piv.final,Pi.final,
                         modBasic,B=B)
# standard errors for the initial probabilities 
options(scipen=999)
round(resboot$sepiv[indC],3)
# standard errors for the transition probabilities from the 1st to the 2nd time occasion
round(resboot$sePi[indC,indC,2],3) 
#

###############################################
source("lmestDecoding.R")
dec <- lmestDecoding(est = mod.final)
# Sort decoded states
Ul2 = dec$Ul
for(u in 1:4) Ul2[dec$Ul==indC[u]] = u
heatmap(Ul2, Rowv = NA, Colv = NA, 
        scale = "none", 
        cexCol = 1.1)
# 
####### Map ################################
####  Map  ####
###############################################
# plot world map by gruop names and not region
require(ggplot2)
require(dplyr)
require(maps)
require(viridis)
require(dichromat)
world_map <- map_data("world")
# region <- factor(world_map$region) # not added

Class2000 <- apply(mod.final$V[,indC,1],1,which.max)
ids<- unique(dt1$id)
ids
classification <- data.frame(group = (ids), 
                             Cluster = Class2000)
row.names(classification) <- NULL
out <- left_join(classification, world_map, by = "group")
colors3 <- rep(c("red", "brown4", "orange", "darkorange1"))
g <- ggplot(out, aes(long, lat, group = Cluster)) +
  geom_polygon(aes( group = group, fill = Cluster)) + 
  theme(legend.position="bottom")
g

#######  ENTROPY ################################
#### Calculate entropy and max entropy  ####
###############################################
#
# 
ent = -sum(mod.final$V*log(pmax(mod.final$V,10^-300))) 
ent
# maximum value
V1 = 0*mod.final$V+1/k.final
entmax = -sum(V1*log(pmax(V1,10^-300))) 
entmax
relenat <- 1-(ent/entmax)
relenat  

#######  Check ################################
#### Check posterior Gaussian distribution ####
###############################################

# posterior probabilities
V0 = out1$mod_HM_now$V
Y <- out1$mod_HM_now$Yimp
require(LMest)
require(otrimle)
# check the size of clusters 
round(out1$mod_HM_now$piv,3)

#### Cluster 1 time 3 #####
# The same should be done for each cluster at each time occasion
# check the distribution of each variable cluster 
ind <- which(t(round(out1$mod_HM_now$V, 0)[, 1, 3]) == 1)
length(ind)
long1 <- matrices2long(Y = Y)
Y_norm_1 <- long1[ind, -c(1, 2)]
# first cluster: measure for each variable
# root mean squared difference between the densities
for(i in 1:4){
  q1<-kerndensmeasure(Y_norm_1[,i])
  print(q1$measure)  
}
names(out1$mod_HM_now$data)


