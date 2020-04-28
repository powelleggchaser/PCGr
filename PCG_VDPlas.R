#Preconditioned Conjugate Gradient (PCG) - Multiple Markers
#Rversion

source("./MME_Iterator_Functions.R")


print_factors = function(n) {
  tmp = rep(0,1)
  #print(paste("The factors of",n,"are:"))
  for(i in 1:n) {
    if((n %% i) == 0) {
     tmp[i] = i
    }
  }
  return(tmp[tmp!="NA"])
}

###########################################################

#### OPTIONS ####

marker_class = "RANDOM" #markers can be treated as fixed or random
genomic_prediction = "NO" #if YES, then validation genotypes read in and GEBVs created for validation set
lamda1=0.1 #shrinkage value applied to model when SNPs fitted as random effects. Not used if marker_class is set to FIXED
relax = 0.8 #relaxation factor to stop divergence of random effect solutions
CONV=1 #set starting convergence value, 1 is suggested
CONV_TOL=0.001 #Set convergence threshold, 10^-9 is suggested

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

#provide empty matrices to store Sum Squares for each Iteration
ssqd_on=as.matrix(rep(0,neq))
ssq_n=as.matrix(rep(0,neq))

#create intital empty vector to store newly generated solutions at each iteration
sol=rep(0,neq)
#create intial empty vector to store previously generated solutions at each iteration
old_sol=rep(0,neq)

C = 1/(t(X)%*%X)
M = diag(neq)
phi = diag(neq)
x = rep(0.1,neq)

r_init = t(X)%*%y - C%*%x
old_r = phi%*%r_init
old_p = rep(0,neq)
old_e = 1

r = rep(0,neq)
p =  rep(0,neq)
e =  0

f50 = print_factors(50)

CONVSampl = rep(NA,10000)
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
    y = (M)%*%old_r
    e = t(y)%*%old_r
    B = e/old_e
    old_e = e
    p[j] = y[j] + (B*old_p[j])
    w = phi%*%(C%*%p)
    
    alpha = e/(t(p)%*%w)
    x = x + as.numeric(alpha)*p[j]
    if (count %in% seq(0,5000,50)){
      r = e - (C%*%x)
    }
    else {r = old_r - as.numeric(alpha)*w[j]}
    
    ssqd_on<-sum(r^2)
    ssq_n<-sum(old_r^2)
    CONV<-ssqd_on/ssq_n
    CONVSampl[j] = CONV
    old_r = r
    old_p = p
    if (CONV<=CONV_TOL) break 
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
