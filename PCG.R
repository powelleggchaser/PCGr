#Preconditioned Conjugate Gradient (PGC) - Multiple Markers
#Rversion

source("./MME_Iterator_Functions.R")

###########################################################

#### OPTIONS ####

marker_class = "RANDOM" #markers can be treated as fixed or random
genomic_prediction = "NO" #if YES, then validation genotypes read in and GEBVs created for validation set
lamda1=0.1 #shrinkage value applied to model when SNPs fitted as random effects. Not used if marker_class is set to FIXED
relax = 0.8 #relaxation factor to stop divergence of random effect solutions
CONV=1 #set starting convergence value, 1 is suggested
CONV_TOL=0.0000000001 #Set convergence threshold, 10^-9 is suggested

###########################################################

### DATA ENTRY ###

###########################################################

#change to directory holding your data files
setwd("./Data_Files/MultipleMarkers")

#read in individuals (Data/No Data Score)
x<-read.csv("X_JH.csv",sep=",",header=F)
x<-as.matrix(x)
class(x)<-"numeric"
ndata=nrow(x)

#read in genotypes
z<-read.table("Genmult_JH.csv",sep=",",header=F)
z<-as.matrix(z)
class(z)<-"numeric"

#read in residuals (required to genertate differing phenotypes)
e<-read.table("E_JH.csv",sep=",",header=F)
e<-as.matrix(e)
class(e)<-"numeric"
average_e<-sum(e)/ndata

#generate phenotypes
allele_effect=read.table("Allele_Effects_JH.csv",sep=",",header=F); allele_effect = as.matrix(allele_effect)
#allele_effect=c(0.5,-1,1,-0.5)
phenotype=rep(0,ndata)
for (i in 1:ndata){
  phenotype[i] = (sum(z[i,]*allele_effect)) + e[i]
}

#############################################################################################################

#calculate allele frequencies & variances
allele_freq_stats=calc_allele_freq_per_locus(z,allele_effect,phenotype)

### MATRIX MANIPULATION ###

#assign phenotypes to y matrix
y=phenotype

#xpy<-as.matrix(t_mult_dupl(x,y))
#zpy<-as.matrix(t_mult_dupl(z,y))

#combine fixed effects and random effects into a single matrix
X=cbind(x,z)

#X<-z
if (marker_class=="FIXED"){
  #Count the number of expected solutions
  neq=ncol(X)
  #Count the number of expected of fixed effect means
  nmeans=ncol(X)
  ######################### SNP EFFECTS = FIXED ##############################
}
if (marker_class=="RANDOM") {
  #Count the number of expected solutions
  neq=ncol(X)
  #Count the number of expected of fixed effect means
  nmeans=ncol(x)
  
  ######################### SNP EFFECTS = RANDOM ##############################
}


r = t(X)%*%y 
C = t(X)%*%X
M = diag(neq)
b = rep(0,neq)
e = r - (C%*%b)
d = solve(M)%*%e

old_r = t(X)%*%y 
old_C = t(X)%*%X
old_M = diag(neq)
old_b = rep(0,neq)
old_e = r - (C%*%b)
old_d = solve(M)%*%e

#provide phenotypes as the starting values for residual updater
#e=as.matrix(y)

#provide empty matrices to store Sum Squares for each Iteration
ssqd_on=as.matrix(rep(0,neq))
ssq_n=as.matrix(rep(0,neq))

#######################################################################################

#################################

####  Preconditioned Conjugate Gradient (PCG)  ####

#################################

#######################################################################################
start = 1
#set counter for iterator
count = 1
#until convergence is reached, iterate process
repeat{
  #for each effect in the model (Fixed + Random)
  for (j in start:neq){
    v = C[j,j]%*%old_d[j]
    w = (t(old_e[j])%*%(solve(M)[j,j]%*%old_e[j])/(t(old_d[j])%*%v))
    b[j] = old_b[j] + w%*%old_d[j]
    e[j] = old_e[j] - w%*%v
    
    B = (t(e[j])%*%v)/(t(old_e[j])%*%(solve(M)[j,j])%*%old_e[j])
    d[j] = v + B%*%old_d[j]
    
    b[j] = b[j] + (relax*(y[j]-C[j,]%*%b))/C[j,j]
    ssqd_on[j]<-sum(e[j]^2)
    ssq_n[j]<-sum(y[j]^2)
    CONV<-ssqd_on[j]/ssq_n[j]
    if (CONV<=CONV_TOL) break 
    old_b[j] = b[j]
    old_d[j] = d[j]
    old_e[j] = e[j]
  }
  count=round(count+1/neq,1)
  print(c("Iteration: ",count))
}


print(paste0("SNP EFFECTS ",marker_class," took ",count," iterations to solve",sep=" "))

sim_effects=allele_effect
plot(sol[-(nmeans)],sim_effects,main="SNP Effect Estimates vs True Values")

ebv<-X%*%sol
gv<-z%*%allele_effect
cor.test(gv,ebv)
plot(gv,ebv,main="Genomic Estimated Breeding Values (GEBVS) vs Genetic Values : Training Set")

######################### TRAINING SET EBV ESTIMATION -  FINISHED ##############################

##########################

### Genomic Prediction ###

##########################

if (genomic_prediction=="YES"){
  #read in validation set individuals (No Phen)
  z_pred=read.csv("Z_Validation.csv",sep=",",header=F);z_pred<-as.matrix(z_pred)
  class(z_pred)<-"numeric"
  ndata=nrow(z_pred)
  effects_estimate<-sol[-(nmeans)]
  
  gebv<-z_pred%*%effects_estimate
  gv = z_pred%*%sim_effects
  cor.test(gv,gebv)
  plot(gv,gebv,main="Genomic Estimated Breeding Values (GEBVS) vs Genetic Values : Validation Set")
  
  ######################### GENOMIC PREDICTION - FINISHED ##############################
}

######################################################################################

### Estimate Variance Components using Bayes C (π=0)

#ss=sum(sol_random^2)+length(e)*Sa
#vara=ss/chi(length(sol_random)+ncol(z))
#ss=sum(e**2)+nue*Se
#vare=ss/chi(nue+ndata) 

######################################################################################
