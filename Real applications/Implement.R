source("Source.R")

################### Data pre processing

#Extract weight matrix, W
hlw <- nb2listw(LO_nb)
hsn <- listw2sn(hlw)
W<- as(hlw, "CsparseMatrix") 

##Extract design matrix, X

hform <- formula(log(price) ~ age + I(age^2) + I(age^3) +
                   + log(lotsize) + rooms + log(TLA) + beds + syear)
hmm0 <- model.matrix(hform, data = house)
X<-hmm0
n<-dim(X)[1]
rownames(W)<-colnames(W)<-1:n
rownames(X)<-1:n


Y<-log(house$price)  # Extract dependent variable, y


## Generate missing values

p_missing <- 0.90     # Proportion of missing values in the response variable (e.g., 0.9 for 90% missing)
p_sample  <- 1 - p_missing  # Proportion of observed (non-missing) values


splitted_sample=splitter(x=X[,-1],y=Y,
                         w=W,p_sample=p_sample)

x.split=splitted_sample$X 
y.split=splitted_sample$Y
w.split=splitted_sample$W
ns=splitted_sample$ns


partineddataset<-partioningdata(x=x.split,y=y.split,w=w.split,no=ns) # For OML

# 1) FML, use complete dataset to estimate parameters


################### Model fitting


tic()
SEM.FML<-ME.full(x=X[,-1],y=Y,
                       w=W,model = 1,Hessian=T,h.full.numerical = F) #SEM
toc()


tic()
SAM.FML<-ME.full(x=X[,-1],y=Y,
                       w=W,model =2,Hessian=T,h.full.numerical = F) #SAM
toc()



round(SEM.FML$Estimates.vec,4)
round(SAM.FML$Estimates.vec,4)

# 2) OML

tic()
SEM.OML<-ME.full(x=partineddataset[[1]],y=partineddataset[[2]],
                       w=partineddataset[[3]],model = 1,Hessian=T,h.full.numerical = F) #SEM
toc()

tic()
SAM.OML<-ME.full(x=partineddataset[[1]],y=partineddataset[[2]],
                       w=partineddataset[[3]],model =2,Hessian=T,h.full.numerical = F) #SAM
toc()


SEM.OML$Estimates.list
SAM.OML$Estimates.list



#3) MML-D


tic()
SEM.MMLD<-ME_DB1(x=x.split,y=y.split,
                           w=w.split,ns=ns,model=1,voo.method ="D", #SEM
                 Hessian = T,h.full.numerical = F)

toc()

tic()
SAM.MMLD<-ME_DB1(x=x.split,y=y.split,
                 w=w.split,ns=ns,model=2,voo.method ="D", #SAM
                 Hessian = T,h.full.numerical = F)

toc()

SEM.MMLD$Estimates.list
SAM.MMLD$Estimates.list

#4 MML-P



tic() 
SEM.MMLP<-ME_B2(x=x.split,y=y.split,w=w.split,no=ns,model=1,Hessian=T, #SEM
                             h.full.numerical = F)
toc()


tic()
SAM.MMLP<-ME_B2(x=x.split,y=y.split,w=w.split,no=ns,model=2,Hessian=T, #SAM
                                        h.full.numerical = F)
toc()


SEM.MMLP$Estimates.list
SAM.MMLP$Estimates.list
