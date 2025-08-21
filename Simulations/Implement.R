source("Source.R")



################ Simulation Settings ################

model <- 2           # Model type: set to 1 for SEM, or 2 for SAM
sqrtn <- 40          # Square root of the number of spatial units (n); e.g., for n = 900, set sqrtn = 30
sqrtn^2              # Total number of spatial units (n = sqrtn^2)

p_missing <- 0.5     # Proportion of missing values in the response variable (e.g., 0.5 for 50% missing)
p_sample  <- 1 - p_missing  # Proportion of observed (non-missing) values

N <- 5              # Number of simulated datasets to generate

# True parameter values (set below) should reflect the underlying data-generating process

beta<-c(1,5) #True values of Betas
rho<-0.8 # True value of rho
sigma2eps<-2 #True value of measurement error variance
sigma2e<-1 #True value of error variance
truepara<-c(beta,rho,sigma2eps,sigma2e)





############### define return lists

OML.estimates.list<-list()
MMLD.estimates.list<-list()
MMLB.estimates.list<-list()



######## RMSE

OML.RMSE<-c()
MMLD.RMSE<-c()
MMLB.RMSE<-c()


#### marginal convergence

# MMLD.conv<-c()
# MMLD.no.iterations<-c()
# 
# MMLB.conv<-c()
# MMLB.no.iterations<-c()
# 
# 
# OML.conv<-c()
# OML.no.iterations<-c()

################### Times

Times<-matrix(rep(0,3*N),ncol=3)

tic() #i=1
for(i in 1:N){
  
  
  
  g<-graph.lattice(c(sqrtn,sqrtn)) 
  g<-connect.neighborhood(g,1) #connect all vertices by one edge.
  wmat<-as_adjacency_matrix(g)
  w<-mat2listw(wmat,style="W") # row normalized w
  w=listw2mat(w)
  w<-as(w,"CsparseMatrix")
  
  nbeta<-2
  p=nbeta-1
  n=ncol(w)
  para<-c(truepara[1:(nbeta+1)],(truepara[length(truepara)])) # With 1-covariates.
  sigma2eps<-truepara[length(truepara)-1]
  
  
  if(model==1){
    sim<-simulateSEM(para,p=p,weightmat=w)
    # model=1
  }else{
    sim<-simulateSAM(para,p=p,weightmat=w)
    # model=2
  }
  
  xsim<-sim$Independent
  ysim<-sim$Dependent
  me=rnorm(n,0,sqrt(sigma2eps))
  
  x=as.matrix(xsim)
  y=ysim+me
  
  splitted_sample=splitter(x=x,y=y,w=w,wmat=wmat,p_sample=p_sample)
  
  
  x=splitted_sample$X 
  y=splitted_sample$Y
  w=splitted_sample$W
  ns=splitted_sample$ns
  
  
  ################################
  
  ########### Fit MLE
  
  Times[i,1]<-system.time(fit.obs.Tom<-ME.full(x=splitted_sample$x_o,y=splitted_sample$y_o,
                          w=splitted_sample$w_oo,
                          model=model,Hessian = T,
                          h.full.numerical = F))[[3]]
  
  
  Times[i,2]<-system.time(fit_marginal_D<-ME_DB1(x=x,y=y,w=w,ns=ns,model=model,voo.method ="D",Hessian=T,
                                                      h.full.numerical = F))[[3]]
  
  Times[i,3]<-system.time(fit_marginal_Bnew<-ME_B2(x=x,y=y,w=w,no=ns,model=model,Hessian=T,
                                                                     h.full.numerical = F))[[3]]
  
  colnames(Times)<-c("OML","MML-D","MML-P")

  
  ###---- collect results---- ####
  
  
  ## OML

  v1<-(diag(solve(fit.obs.Tom$Hess.mat)))
  estimates.vec<-fit.obs.Tom$Estimates.vec
  ci1<-createCI_1(variances=v1,estimate=fit.obs.Tom$Estimates.vec,truepara=truepara)
  OML.estimates.list[[i]]<-ci1
  
  ## MML-B
  
  v2<-(diag(solve(fit_marginal_Bnew$Hessian)))
  ci2<-createCI_1(variances=v2,estimate=fit_marginal_Bnew$Estimates.vec,truepara=truepara)
  MMLB.estimates.list[[i]]<-ci2
  
  ### MML-D

  vmarginal<-(diag(solve(fit_marginal_D$Hessian)))
  ci_marginal<-createCI_1(variances=vmarginal,estimate=fit_marginal_D$Estimates.vec,truepara=truepara)
  MMLD.estimates.list[[i]]<-ci_marginal
  # MMLD.conv[i]<-fit_marginal_D$optim.result$convergence

  
  ### RMSE calculation
  nom.p<-length(truepara)
  MMLD.RMSE[i]<-sqrt((1/nom.p)*sum((fit_marginal_D$Estimates.vec-truepara)^2))
  OML.RMSE[i]<-sqrt((1/nom.p)*sum((estimates.vec-truepara)^2))
  MMLB.RMSE[i]<-sqrt((1/nom.p)*sum((fit_marginal_Bnew$Estimates.vec-truepara)^2))
  
  print(i)
  
}
toc()

############# Summarize results across datasets.

# Accuracy


###############################

#Mean estimates, and MSE across simulated datasets

(Reduce(`+`, OML.estimates.list)/N )[,1:3] #OML
(Reduce(`+`, MMLD.estimates.list)/N)[,1:3] #MMLD
(Reduce(`+`, MMLB.estimates.list)/N)[,1:3] #MMLB



# RMSE

mean(OML.RMSE)
mean(MMLD.RMSE)
mean(MMLB.RMSE)


##### Remove results for datasets with numerical issues in variance calculation

# Find estimates (induces) with all NA in coverage.


vec.NA.or.valid<-unlist(lapply(OML.estimates.list, function(est.mat) sum(est.mat[,7])))
vec.NA.or.valid.OLM <- !is.na(vec.NA.or.valid)
sum(vec.NA.or.valid.OLM)

vec.NA.or.valid<-unlist(lapply(MMLD.estimates.list, function(est.mat) sum(est.mat[,7])))
vec.NA.or.valid.marginalD <- !is.na(vec.NA.or.valid)
sum(vec.NA.or.valid.marginalD)


vec.NA.or.valid<-unlist(lapply(MMLB.estimates.list, function(est.mat) sum(est.mat[,7])))
vec.NA.or.valid.marginalB <- !is.na(vec.NA.or.valid)
sum(vec.NA.or.valid.marginalB)



final.valid.indeces<- vec.NA.or.valid.OLM & vec.NA.or.valid.marginalD & vec.NA.or.valid.marginalB

## Summary of data sets with var is not NA

(Reduce(`+`, OML.estimates.list[final.valid.indeces]) / sum(final.valid.indeces))
(Reduce(`+`, MMLD.estimates.list[final.valid.indeces]) / sum(final.valid.indeces))
(Reduce(`+`, MMLD.estimates.list[final.valid.indeces]) / sum(final.valid.indeces))




# Average computation times

apply(Times, 2, mean)
boxplot(Times)




