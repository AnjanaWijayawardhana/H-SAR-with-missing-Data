library(numDeriv)
library(Matrix)
library(spatialreg)
library(spData)
library(spdep)
library(tictoc)
library(igraph)





#################### MLE algorithms #####################

ME.full<-function(x,y,w,model,Hessian,
                  h.full.numerical){
  
  t1<-Sys.time()
  w<-as(w,"CsparseMatrix")
  n=ncol(w)
  
  w=Matrix(w) ## Matrix
  I=Diagonal(n) #Matrix
  
  wpluswt=(w+t(w))
  wwt=w%*%t(w)
  
  x=cbind(rep(1,n),x)
  
  WX<-w%*%x
  WY<-w%*%y  
  
  # para=c(0.1,0.2)
  log.like.rho.theta<-function(para,w,x,y,WX,WY,I,model,wpluswt,wwt){
    
    n=nrow(w)
    
    rho=para[1]
    theta=exp(para[2])
    
    AY<- y-rho*WY
    
    
    if(model==1){
      AX<-x-rho*WX
    }else{
      AX<-x
    }
    
    
    K<-I-rho*wpluswt+rho^2*wwt  #A*A' =(I-rho*W)*(I-rho*W)' I -rho*(W+W') +rho^2*(W*W')
    K<-forceSymmetric(K)
    
    
    
    chol.K <- Cholesky(K, super = TRUE, Imult = 0)
    
    ldet.K<- 2*determinant(chol.K)$modulus #compare log(det(K))
    
    
    
    
    chol.K.theta<- update(chol.K, K, mult=theta)
    ldet.K.theta<-2*determinant(chol.K.theta)$modulus
    logdetV<- ldet.K.theta - ldet.K  #compare with log(det(I+theta*solve((I-rho*t(W))%*%(I-rho*W))))
    
    c<- solve(chol.K.theta,AY,system = "A")
    C<- solve(chol.K.theta,AX,system = "A")
    
    
    CtC<- t(AX)%*%C
    Ctc<- t(AX)%*%c  # t(AX)%*%solve(K+theta*I)%*%AY # compare with Ctb
    
    ctc<- t(AY)%*%c  # t(AY)%*%solve(K+theta*I)%*%AY # compare with btb
    beta<-solve(CtC)%*%Ctc
    rtMr<-as.numeric(ctc - t(Ctc)%*%beta)  
    
    omega<-rtMr/n
    
    # L<- c(0.5*n*log(2*pi) + 0.5*n1*log(omega)   + 0.5*logdetV + 0.5*rtMr/omega + add)
    
    L<--(-0.5*n*log(2*pi)-0.5*n*log(omega)-0.5*logdetV-n*0.5) #
    L<-as.numeric(L)
    
  }
  
  
  
  
  para=c(0.1,log(0.1))
  optimresult=optim(par=para,fn=log.like.rho.theta,w=w,x=x,y=y,WX=WX,WY=WY,I=I,model=model,wpluswt=wpluswt,wwt=wwt,
                    method = "L-BFGS-B",lower = c(-0.99999,-10), upper = c(0.99999,10),
                    control = list(trace = F))
  
  
  rho_theta=optimresult$par
  count=optimresult$counts
  
  rho=rho_theta[1]
  theta=exp(rho_theta[2])
  # rho=0.1
  # theta=0.1
  ####
  
  AY<- y-rho*WY
  if(model==1){
    AX<-x-rho*WX
  }else{
    AX<-x
  }
  
  
  K<-I-rho*wpluswt+rho^2*wwt  #A*A' =(I-rho*W)*(I-rho*W)' I -rho*(W+W') +rho^2*(W*W')
  K<-forceSymmetric(K)
  
  
  
  chol.K <- Cholesky(K, super = TRUE, Imult = 0)
  
  ldet.K<- 2*determinant(chol.K)$modulus #compare log(det(K))
  
  chol.K.theta<- update(chol.K, K, mult=theta)
  ldet.K.theta<-2*determinant(chol.K.theta)$modulus
  logdetV<- ldet.K.theta - ldet.K  #compare with log(det(I+theta*solve((I-rho*t(W))%*%(I-rho*W))))
  
  c<- solve(chol.K.theta,AY,system = "A")
  C<- solve(chol.K.theta,AX,system = "A")
  
  CtC<- t(AX)%*%C
  Ctc<- t(AX)%*%c  # t(AX)%*%solve(K+theta*I)%*%AY # compare with Ctb
  
  
  
  ctc<- t(AY)%*%c  # t(AY)%*%solve(K+theta*I)%*%AY # compare with btb
  beta<-solve(CtC)%*%Ctc
  rtMr<-as.numeric(ctc - t(Ctc)%*%beta)  
  
  omega<-rtMr/n
  
  ####
  
  sigma2eps<-omega
  sigma2y<-theta*omega
  
  estimates.vec<-c(as.vector(beta),as.numeric(rho),
                   as.numeric(sigma2eps),as.numeric(sigma2y))
  
  estimates.list<-list(beta=(beta),rho=as.numeric(rho),
                       sigma2eps=as.numeric(sigma2eps),sigma2y=as.numeric(sigma2y))
  
  t2<-Sys.time()
  Est.time<-t2-t1
  
  log_like_full<-function(para,X,y,I,w,model){
    
    x_till<-X
    beta<-para[1:ncol(x_till)]
    rho<-para[ncol(x_till)+1]
    sigma2eps<-para[ncol(x_till)+2]
    sigma2y=para[ncol(x_till)+3]
    
    theta<-sigma2y/sigma2eps
    
    if(theta<0){
      return(1e6) #to avoid NA in log(sigma2y), we lilit that theta should be positive
    }
    
    A<-I-rho*w
    
    
    AtA=t(A)%*%A
    v=(I+theta*solve(AtA))
    
    V=forceSymmetric(v) # make symetric
    # det(V_oo)
    if(model==1){
      r<-y-(x_till%*%beta)
    }else{
      r<-y-solve(A,x_till%*%beta)
    }
    
    
    logdet<-determinant((V),logarithm=T)$modulus[[1]]
    
    L<--0.5*n*log(2*pi)-0.5*n*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(r)%*%solve(V,r)
    return(as.numeric(-L))
    
  } # My log-like full variables
  
  log_like_rho.sigmas<-function(para,beta,X,y,I,w,model){
    
    x_till<-X
    
    rho<-para[1]
    sigma2eps<-para[2]
    sigma2y=para[3]
    
    theta<-sigma2y/sigma2eps
    
    if(theta<0){
      return(1e6) #to avoid NA in log(sigma2y), we lilit that theta should be positive
    }
    
    A<-I-rho*w
    
    
    AtA=t(A)%*%A
    v=(I+theta*solve(AtA))
    
    V=forceSymmetric(v) # make symetric
    # det(V_oo)
    if(model==1){
      r<-y-(x_till%*%beta)
    }else{
      r<-y-solve(A,x_till%*%beta)
    }
    
    
    logdet<-determinant((V),logarithm=T)$modulus[[1]]
    
    L<--0.5*n*log(2*pi)-0.5*n*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(r)%*%solve(V,r)
    return(as.numeric(-L))
    
  } # My log-like rho.sigmas func
  
  
  log_like_rho.sigmas_eff<-function(para,beta,X,y,I,w,model){
    
    x_till<-X
    rho<-para[1]
    sigma2eps<-para[2]
    sigma2y=para[3]
    
    theta<-sigma2y/sigma2eps
    
    if(theta<0){
      return(1e6) #to avoid NA in log(sigma2y), we lilit that theta should be positive
    }
    
    A<-I-rho*w
    
    
    AtA=t(A)%*%A
    
    v=(I+theta*solve(Cholesky(forceSymmetric(AtA)),I))
    # ?solve
    V=forceSymmetric(v) # make symetric
    # det(V_oo)
    if(model==1){
      r<-y-(x_till%*%beta)
    }else{
      r<-y-solve(A,x_till%*%beta)
    }
    
    
    logdet<-determinant((V),logarithm=T)$modulus[[1]]
    
    L<--0.5*n*log(2*pi)-0.5*n*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(r)%*%solve(V,r)
    return(as.numeric(-L))
    
  } # My log-like rho.sigmas func_eff
  
  
  cal.info.beta<-function(x,y,I,w,wpluswt,wwt,rho,sigma2eps,sigma2y,model){
    
    theta<-sigma2y/sigma2eps
    K<-I-rho*wpluswt+rho^2*wwt
    M<-solve(I+theta*solve(K))
    if(model==1){
      return((t(x)%*%M%*%x)/sigma2eps)
    }else{
      A<-I-rho*w
      sol_Ax<-solve(A,x)
      dr_by_dbeta<-0-sol_Ax
      return(-(1/sigma2eps)*t(sol_Ax)%*%M%*%dr_by_dbeta)
    }
    
  } # cal info beta
  
  cal.info.beta_eff<-function(x,y,I,w,wpluswt,wwt,rho,sigma2eps,sigma2y,model){
    
    theta<-sigma2y/sigma2eps
    K<-I-rho*wpluswt+rho^2*wwt
    if(model==1){
      M<-solve(I+theta*solve(K))
      return((t(x)%*%M%*%x)/sigma2eps)
    }else{
      return((1/sigma2eps)*t(x)%*%solve(K+theta*I,x))
    }
    
  } # cal info beta eff
  
  
  
  t3<-Sys.time()
  if(Hessian==T){
    
    if(h.full.numerical){
      Hess.mat<-hessian(func=log_like_full,x=estimates.vec,I=I,y=y,w=w,X=x,model=model)
    }else{
      info.beta<-cal.info.beta_eff(x=x,y=y,I=I,w=w,wpluswt=wpluswt,wwt=wwt,
                                   rho=rho,sigma2eps=sigma2eps,sigma2y=sigma2y,model=model)
      
      hesmat_rho.sigmas_eff<-hessian(func=log_like_rho.sigmas_eff,x=estimates.vec[(ncol(x)+1):length(estimates.vec)],
                                     beta=as.vector(beta),I=I,y=y,w=w,X=x,model=model)
      
      hesmatrhosigmas_3<-0
      hesmatrhosigmas_4<-0
      # hesmatrhosigmas<-0
      p<-length(as.vector(beta))
      Hess.mat<-matrix(rep(0,length(estimates.vec)*length(estimates.vec)),ncol=length(estimates.vec))
      Hess.mat[1:p,1:p]<-as.matrix(info.beta)
      Hess.mat[(1+p):length(estimates.vec),(1+p):length(estimates.vec)]<-
        hesmat_rho.sigmas_eff
    }
    
    
  }else{
    Hess.mat<-0
    hesmat_rho.sigmas_eff<-0
    info.beta<-0
    hesmat_rho.sigmas_eff<-0
  }
  t4<-Sys.time()
  Hess.time<-t4-t3
  
  
  Times<-c(Est.time,Hess.time)
  names(Times)<-c("Est.time","Hess.time")
  
  
  return(list(Estimates.vec=estimates.vec,
              Estimates.list=estimates.list,Times=Times,no.of.function.eval=count,
              Hess.mat=Hess.mat,optim.result=optimresult))
  
  
} # MLE of Suesse (2018): MLE for SAR with measurement error. 
                                        #This algorithm is used to OML using CCA.


ME_DB1<-function(x,y,w,ns,model,voo.method,
               Hessian,h.full.numerical){ #include
  
  w<-as(w,"CsparseMatrix")
  t1<-Sys.time()
  y<-as.matrix(y,ncol=1)
  w=Matrix(w)  #Matrix
  n=nrow(w)
  I=Diagonal(n) #Matrix
  nu=n-ns
  x<-cbind(1,x)
  
  
  
  if(voo.method=="B"){
    b1s_11=Diagonal(ns)
    b1s_12=Matrix(0,ns,nu)
    b1s=cbind(b1s_11,b1s_12)
    
    wpluswt=(w+t(w))
    wwt=w%*%t(w)
    wtw=t(w)%*%w
    log_like<-function(par,x,y,w,ns,I,b1s,model,wpluswt,wwt){ # I -Matrix, b1s-Matrix,w, Matrix
      
      rho=par[1]
      theta=exp(par[2])
      n=nrow(w)
      
      A=I-rho*w
      AAt=(I-rho*(wpluswt)+rho*rho*wwt)
      #vss
      b=solve(t(A),t(b1s))
      vss=forceSymmetric(t(b)%*%(AAt+theta*I)%*%b)
      
      
      if(model==1){
        x_till=x
      }else{
        x_till=solve(A,x)
      }
      
      x_till_s=x_till[(1:ns),]
      y_s=as.matrix(y[(1:ns),])
      
      beta_hat=solve(t(x_till_s)%*%solve(vss,x_till_s))%*%t(x_till_s)%*%solve(vss,y_s)
      #beta_hat=c(1,1,2,2)
      #############################################
      # r=y-x_till%*%beta_hat
      # rs=r[1:ns,]
      
      rs<-y[1:ns,]-x_till_s%*%beta_hat
      ###############################################
      #b=solve(t(A),t(b1s))
      #vss=t(b)%*%(A%*%t(A)+theta*I)%*%b
      rt_vssin_r=t(rs)%*%solve(vss,rs)
      
      ##########log det part##########
      vss=forceSymmetric(vss) #To make vss symmetric
      
      vss_sparse=vss
      
      log_det=-sum(log(diag(chol(vss_sparse))))
      
      l=-0.5*ns*log(2*pi)-0.5*ns*log(rt_vssin_r/ns)+log_det-0.5*ns
      l=as.numeric(l)
      return(-l)
    }
    
    para=c(0.1,log(0.1))
    
    optimresult=optim(par=para,fn=log_like,x=x,y=y,w=w,ns=ns,I=I,b1s=b1s,model=model,wpluswt=wpluswt,wwt=wwt,method = "L-BFGS-B",
                      lower = c(-0.99999,-10), upper = c(0.99999,10),hessian = F)
    
    
    LL.val<-optimresult$value
    rho_theta=optimresult$par
    count=optimresult$counts
    # optimresult$hessian
    
    # estimates
    rho=rho_theta[1]
    theta=exp(rho_theta[2])
    rho_hat=rho
    theta_hat=theta
    
    A=I-rho*w
    
    #vss
    b=solve(t(A),t(b1s))
    vss=forceSymmetric(t(b)%*%(A%*%t(A)+theta*I)%*%b)
    
  }else{
    wtplusw=(w+t(w))
    wtw=t(w)%*%(w)
    log_like<-function(par,x,y,w,ns,I,model,wtplusw,wtw){ # I -Matrix, b1s-Matrix,w, Matrix
      
      rho=par[1]
      theta=exp(par[2])
      n=nrow(w)
      
      A=I-rho*w
      AtA=(I-rho*(wtplusw)+rho*rho*wtw)
      L_ATA<-Cholesky(AtA)
      vss<-(I+theta*(solve(L_ATA,I)))[1:ns,1:ns]
      vss=forceSymmetric(vss)
      
      if(model==1){
        x_till=x
      }else{
        x_till=solve(A,x)
      }
      
      x_till_s=x_till[(1:ns),]
      y_s=as.matrix(y[(1:ns),])
      
      beta_hat=solve(t(x_till_s)%*%solve(vss,x_till_s))%*%t(x_till_s)%*%solve(vss,y_s)
      #beta_hat=c(1,1,2,2)
      
      rs<-y[1:ns,]-x_till_s%*%beta_hat
      ###############################################
      #b=solve(t(A),t(b1s))
      #vss=t(b)%*%(A%*%t(A)+theta*I)%*%b
      rt_vssin_r=t(rs)%*%solve(vss,rs)
      
      ##########log det part######
      vss_sparse=vss
      log_det=-sum(log(diag(chol(vss_sparse))))
      
      l=-0.5*ns*log(2*pi)-0.5*ns*log(rt_vssin_r/ns)+log_det-0.5*ns
      l=as.numeric(l)
      return(-l)
    }
    para=c(0.1,log(0.1))
    optimresult=optim(par=para,fn=log_like,x=x,y=y,w=w,ns=ns,I=I,
                      model=model,wtplusw=wtplusw,wtw=wtw,method = "L-BFGS-B",
                      lower = c(-0.99999,-10), upper = c(0.99999,10),hessian = F)
    
    
    LL.val<-optimresult$value
    rho_theta=optimresult$par
    count=optimresult$counts
    # optimresult$hessian
    
    # estimates
    rho=rho_theta[1]
    theta=exp(rho_theta[2])
    rho_hat=rho
    theta_hat=theta
    
    A=I-rho*w
    AtA=(I-rho*(wtplusw)+rho*rho*wtw)
    L_ATA<-Cholesky(AtA)
    vss<-(I+theta*(solve(L_ATA,I)))[1:ns,1:ns]
    vss=forceSymmetric(vss)
    
    
  }
  
  
  #betahat and r
  
  
  if(model==1){
    x_till=x
  }else{
    x_till=solve(A,x)
  }
  
  x_till_s=x_till[(1:ns),]
  y_s=as.matrix(y[(1:ns),])
  
  beta_hat=solve(t(x_till_s)%*%solve(vss,x_till_s))%*%t(x_till_s)%*%solve(vss,y_s)
  
  rs<-y[1:ns,]-x_till_s%*%beta_hat
  
  rt_vssin_r=t(rs)%*%solve(vss,rs)
  omega_hat=rt_vssin_r/ns
  
  
  log_like_full<-function(para,X,y,I,w,no){
    # no=ns
    # para<-c(1,2,0.1,0.3,0.4)
    x<-X
    xo<-x[1:no,]
    yo<-y[1:no,]
    
    beta<-para[1:ncol(xo)]
    rho<-para[ncol(xo)+1]
    sigma2eps<-para[ncol(xo)+2]
    sigma2y=para[ncol(xo)+3]
    
    theta<-sigma2y/sigma2eps
    
    
    A<-I-rho*w
    AtA=t(A)%*%A
    v=(I+theta*solve(AtA))
    V=forceSymmetric(v) # make symetric
    V_oo=V[1:no,1:no]
    # det(V_oo)
    ro<-yo-xo%*%beta
    
    logdet<-determinant((V_oo),logarithm=T)$modulus[[1]]
    
    L<--0.5*no*log(2*pi)-0.5*no*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(ro)%*%solve(V_oo,ro)
    return(as.numeric(-L))
    
  } # My marginal log-like full variables
  
  
  log_like_full<-function(para,X,y,I,w,no,model){
    # no=ns
    # para<-c(1,2,0.1,0.3,0.4)
    x<-X
    xo<-x[1:no,]
    yo<-y[1:no,]
    
    beta<-para[1:ncol(xo)]
    rho<-para[ncol(xo)+1]
    sigma2eps<-para[ncol(xo)+2]
    sigma2y=para[ncol(xo)+3]
    
    theta<-sigma2y/sigma2eps
    
    
    A<-I-rho*w
    AtA=t(A)%*%A
    v=(I+theta*solve(AtA))
    V=forceSymmetric(v) # make symetric
    V_oo=V[1:no,1:no]
    # det(V_oo)
    if(model==1){
      ro<-yo-xo%*%beta
    }else{
      ro<-yo-(solve(A,x)[(1:no),])%*%beta
    }
    
    
    logdet<-determinant((V_oo),logarithm=T)$modulus[[1]]
    
    L<--0.5*no*log(2*pi)-0.5*no*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(ro)%*%solve(V_oo,ro)
    return(as.numeric(-L))
    
  } # My marginal log-like full variables
  
  log_like_B<-function(para,beta,X,y,I,w,no,
                       b1s,wpluswt,wwt){
    # no=ns
    # para<-c(1,2,0.1,0.3,0.4)
    x<-X
    xo<-x[1:no,]
    yo<-y[1:no,]
    
    
    rho<-para[1]
    sigma2eps<-para[2]
    sigma2y=para[3]
    
    theta<-sigma2y/sigma2eps
    
    if(theta<0){
      return(1e6) #to avoid NA in log(sigma2y), we lilit that theta should be positive
    }
    
    A=I-rho*w
    AAt=(I-rho*(wpluswt)+rho*rho*wwt)
    #vss
    b=solve(t(A),t(b1s))
    V_oo=forceSymmetric(t(b)%*%(AAt+theta*I)%*%b)
    
    ro<-yo-xo%*%beta
    
    logdet<-determinant((V_oo),logarithm=T)$modulus[[1]]
    
    L<--0.5*no*log(2*pi)-0.5*no*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(ro)%*%solve(V_oo,ro)
    return(as.numeric(-L))
    
  } # correct efficient negative log likelihood of rho,sigmas- B method
  
  
  log_like_D<-function(para,beta,X,y,I,w,no,wpluswt,wtw){
    # no=ns
    # para<-c(1,2,0.1,0.3,0.4)
    x<-X
    xo<-x[1:no,]
    yo<-y[1:no,]
    
    
    rho<-para[1]
    sigma2eps<-para[2]
    sigma2y=para[3]
    
    theta<-sigma2y/sigma2eps
    
    
    A<-I-rho*w
    AtA=(I-rho*(wpluswt)+rho*rho*wtw)
    
    v=(I+theta*solve(Cholesky(forceSymmetric(AtA)),I))
    V=forceSymmetric(v) # make symetric
    V_oo=V[1:no,1:no]
    # det(V_oo)
    
    # beta=solve(t(xo)%*%solve(V_oo,xo))%*%t(xo)%*%solve(V_oo,yo)
    
    ro<-yo-xo%*%beta
    
    logdet<-determinant((V_oo),logarithm=T)$modulus[[1]]
    
    L<--0.5*no*log(2*pi)-0.5*no*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(ro)%*%solve(V_oo,ro)
    return(as.numeric(-L))
    
  } # correct inefficient negative log likelihood of rho,sigmas- D method
  
  
  sigma2eps<-omega_hat
  sigma2y<-theta_hat*omega_hat
  
  estimates<-c(as.vector(beta_hat),as.numeric(rho_hat),
               as.numeric(sigma2eps),as.numeric(sigma2y))
  
  t2<-Sys.time()
  if(Hessian){
    
    if(h.full.numerical){
      hesmat_full<-numDeriv::hessian(func=log_like_full,
                                     x=estimates,
                                     I=I,X=x,y=y,w=w,no=ns,
                                     model=model)
      # diag(solve(hesmat_full))
    }else{
      if(voo.method!="B"){
        # print("Direct method is used when calculating Hessian.....")
        wpluswt=(t(w)+w)
        wtw=t(w)%*%(w)
        Hess.rho.sigmas<-numDeriv::hessian(func=log_like_D,beta=as.matrix(beta_hat,ncol=1),
                                           x=estimates[(ncol(x)+1):length(estimates)],
                                           I=I,X=x_till,y=y,wpluswt=wpluswt,
                                           wtw=wtw,w=w,no=ns,method.args=list(eps=1e-4,
                                                                              d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
        
        
      }else{
        # print("B method is used when calculating Hessian.....")
        wpluswt=(w+t(w))
        wwt=w%*%t(w)
        Hess.rho.sigmas<-numDeriv::hessian(func=log_like_B,beta=as.matrix(beta_hat,ncol=1),
                                           x=estimates[(ncol(x)+1):length(estimates)],
                                           I=I,X=x_till,
                                           b1s=b1s,y=y,wpluswt=wpluswt,
                                           wwt=wwt,w=w,no=ns,method.args=list(eps=1e-4,
                                                                              d=0.0001,
                                                                              zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
        
      }
      
      nofbetas<-ncol(x_till_s)
      info<-matrix(rep(0,(ncol(x)+3)^2),ncol=ncol(x)+3)
      
      exactinfo.beta<-as((t(x_till_s)%*%solve(vss,x_till_s))*as.numeric((1/sigma2eps)),
                         "matrix")
      info[1:nofbetas,1:nofbetas]<-exactinfo.beta
      info[(nofbetas+1):nrow(info),(nofbetas+1):ncol(info)]<-Hess.rho.sigmas
      hesmat_full<-info
    }
    
    
    
    
  }else{
    hesmat_full<-0
    hess_2<-0
  }
  t3<-Sys.time()
  Hess.time<-t3-t2
  Est.time<-t2-t1
  Times<-c(Est.time,Hess.time)
  names(Times)<-c("Est.time","Hess.time")
  
  estimates.list=list(beta=beta_hat,rho=rho_hat,sigma2eps=sigma2eps,sigma2y=sigma2y)
  
  # return(list(Estimates=l,Hessian=hesmat_full,Times=Times))
  info.beta.rho<-(t(solve(A,x)[1:ns,])%*%solve(vss)%*%((solve(A,w)%*%solve(A,x))[1:ns,])%*%beta_hat)/as.numeric(sigma2eps)
  return(list(Estimates.vec=estimates,Hessian=hesmat_full,
              Estimates.list=estimates.list,Times=Times,optim.result=optimresult,info.beta.rho=info.beta.rho))
  
} # MML-D algorithm,



ME_B2<-function(x,y,w,no,model,  
                           Hessian,h.full.numerical){ #include
  
  w<-as(w,"CsparseMatrix")
  t1<-Sys.time()
  y<-as.matrix(y,ncol=1)
  w=Matrix(w)  #Matrix
  n=nrow(w)
  I=Diagonal(n) #Matrix
  ns=no
  nu=n-ns
  x<-cbind(1,x)
  
  
  b1s_11=Diagonal(ns)
  b1s_12=Matrix(0,ns,nu)
  b1s=cbind(b1s_11,b1s_12)
  Bo=b1s
  no=ns
  yo=as.matrix(y[(1:no),])
  yotyo<-t(yo)%*%yo
  yotBo<-t(yo)%*%Bo
  BotBo=t(Bo)%*%Bo
  
  wpluswt=(w+t(w))
  wwt=w%*%t(w)
  
  wtw=t(w)%*%w
  
  log_like2<-function(par,x,yo,yotyo,yotBo,BotBo,w,no,I,Bo,model,wpluswt,wtw){ # I -Matrix, b1s-Matrix,w, Matrix
    
    rho=par[1]
    theta=exp(par[2])
    n=nrow(w)
    
    A=I-rho*w
    AtA=(I-rho*(wpluswt)+rho*rho*wtw)
    #vss
    
    
    if(model==1){
      x_till=x
    }else{
      x_till=solve(A,x)
    }
    
    x_till_o=x_till[(1:no),]
    
    
    # chol.AtAplus=Cholesky(AtA+theta*t(Bo)%*%(Bo))
    # beta_hat_1=solve(t(x_till_o)%*%x_till_o-theta*t(x_till_o)%*%Bo%*%solve(chol.AtAplus,t(Bo)%*%x_till_o))
    # beta_hat_2=t(x_till_o)%*%yo-theta*t(x_till_o)%*%Bo%*%solve(chol.AtAplus,t(Bo)%*%yo)
    # beta_hat_1%*%beta_hat_2
    
    
    x_til_otx_till_o=t(x_till_o)%*%x_till_o
    x_til_otyo<-t(x_till_o)%*%yo
    
    chol.AtAplus=Cholesky(AtA+theta*BotBo)
    beta_hat_1=solve(x_til_otx_till_o-theta*t(x_till_o)%*%Bo%*%solve(chol.AtAplus,t(Bo)%*%x_till_o))
    beta_hat_2=x_til_otyo-theta*t(x_till_o)%*%Bo%*%solve(chol.AtAplus,t(yotBo))
    beta_hat=beta_hat_1%*%beta_hat_2
    
    
    
    # ro<-yo-x_till_o%*%beta_hat
    # Bot.ro<-t(Bo)%*%ro
    # rotvooinvro=t(ro)%*%ro-theta*t(Bot.ro)%*%solve(chol.AtAplus,Bot.ro)
    
    
    rotvooinvro=yotyo-2*t(beta_hat)%*%x_til_otyo+t(beta_hat)%*%x_til_otx_till_o%*%beta_hat-
      theta*(yotBo-t(beta_hat)%*%t(x_till_o)%*%Bo)%*%solve(chol.AtAplus,(t(Bo)%*%yo-t(Bo)%*%x_till_o%*%beta_hat))
    
    # sum(theta*t(Bot.ro))
    # sum(theta*(yotBo-t(beta_hat)%*%t(x_till_o)%*%Bo))
    
    ##########log det part##########
    
    # method 1
    # tic()
    # AtA_L<-Cholesky(AtA)
    # BoBot=Bo%*%t(Bo)
    # voo1=BoBot+theta*Bo%*%solve(AtA_L,t(Bo))
    # log_det=-sum(log(diag(chol(voo1))))
    # toc()
    
    # log(det(AtA+theta*t(Bo)%*%Bo))
    
    # log(det(AtA))+log(det(I+theta*t(Bo)%*%Bo))
    # determinant(AtA+theta*t(Bo)%*%Bo)$modulus
    
    # log_det=-0.5*(log(det(Diagonal(no)))-log(det(AtA))+log(det(AtA+theta*t(Bo)%*%Bo)))
    # log_det1=-0.5*(-determinant(AtA)$modulus+determinant(AtA+theta*t(Bo)%*%Bo)$modulus)
    
    # method 2
    # tic()
    AtA_L<-Cholesky(AtA)
    log_det1=-0.5*(-determinant(AtA_L)$modulus*2+determinant(chol.AtAplus)$modulus*2) # cal det using cholfac
    log_det1=as.numeric(log_det1)
    # toc()
    
    
    
    
    # chol.K.theta<- update(chol.K, K, mult=theta)
    # chol.AtAplus2=update(AtA_L,theta*t(Bo)%*%Bo)
    # class(AtA_L)
    # determinant(chol.AtAplus2)$modulus*2
    
    l=-0.5*ns*log(2*pi)-0.5*no*log(rotvooinvro/no)+log_det1-0.5*no
    l=as.numeric(l)
    return(-l)
  }
  
  
  
  para=c(0.1,log(0.1))
  
  optimresult=optim(par=para,fn=log_like2,x=x,y=yo,yotyo=yotyo,yotBo=yotBo,BotBo=BotBo,w=w,no=no,
                    I=I,Bo=Bo,model=model,wpluswt=wpluswt,wtw=wtw,method = "L-BFGS-B",
                    lower = c(-0.99999,-10), upper = c(0.99999,10),hessian = F)
  LL.val<-optimresult$value
  rho_theta=optimresult$par
  count=optimresult$counts
  # optimresult$hessian
  
  # estimates
  rho=rho_theta[1]
  theta=exp(rho_theta[2])
  rho_hat=rho
  theta_hat=theta
  
  A=I-rho*w
  AtA=(I-rho*(wpluswt)+rho*rho*wtw)
  #vss
  
  #vss
  # b=solve(t(A),t(b1s))
  # vss=forceSymmetric(t(b)%*%(A%*%t(A)+theta*I)%*%b)
  
  #betahat and r
  
  
  if(model==1){
    x_till=x
  }else{
    x_till=solve(A,x)
  }
  
  x_till_o=x_till[(1:no),]
  y_o=as.matrix(y[(1:no),])
  
  chol.AtAplus=Cholesky(AtA+theta*t(Bo)%*%(Bo))
  
  beta_1_inner=t(x_till_o)%*%x_till_o-theta*t(x_till_o)%*%Bo%*%solve(chol.AtAplus,t(Bo)%*%x_till_o)
  beta_hat_1=solve(beta_1_inner)
  beta_hat_2=t(x_till_o)%*%y_o-theta*t(x_till_o)%*%Bo%*%solve(chol.AtAplus,t(Bo)%*%y_o)
  beta_hat=beta_hat_1%*%beta_hat_2
  
  
  
  ro<-y_o-x_till_o%*%beta_hat
  Bot.ro<-t(Bo)%*%ro
  solve.chol.AtAplus.Bot.ro<-solve(chol.AtAplus,Bot.ro)
  rotvooinvro=t(ro)%*%ro-theta*t(Bot.ro)%*%solve.chol.AtAplus.Bot.ro
  omega_hat=rotvooinvro/no
  
  
  
  log_like_B<-function(para,beta,X,yo,yotyo,yotBo,BotBo,I,w,no,
                       Bo,wpluswt,wtw){
    # no=ns
    # para<-c(1,2,0.1,0.3,0.4)
    x<-X
    xo<-x[1:no,]
    # yo<-y[1:no,]
    
    
    rho<-para[1]
    sigma2eps<-para[2]
    sigma2y=para[3]
    
    theta<-sigma2y/sigma2eps
    
    if(theta<0){
      return(1e6) #to avoid NA in log(sigma2y), we lilit that theta should be positive
    }
    
    A=I-rho*w
    AtA=(I-rho*(wpluswt)+rho*rho*wtw)
    #vss
    # b=solve(t(A),t(b1s))
    # V_oo=forceSymmetric(t(b)%*%(AAt+theta*I)%*%b)
    # ro<-yo-xo%*%beta
    # logdet<-determinant((V_oo),logarithm=T)$modulus[[1]]
    # 
    # 
    x_til_otyo<-t(x_till_o)%*%yo
    x_til_otx_till_o=t(x_till_o)%*%x_till_o
    chol.AtAplus=Cholesky(AtA+theta*t(Bo)%*%(Bo))
    rotvooinvro=yotyo-2*t(beta)%*%x_til_otyo+t(beta)%*%x_til_otx_till_o%*%beta-
      theta*(yotBo-t(beta)%*%t(x_till_o)%*%Bo)%*%solve(chol.AtAplus,(t(Bo)%*%yo-t(Bo)%*%x_till_o%*%beta))
    
    AtA_L<-Cholesky(AtA)
    log_det1=-0.5*(-determinant(AtA_L)$modulus*2+determinant(chol.AtAplus)$modulus*2) # cal det using cholfac
    log_det1=as.numeric(log_det1)
    
    
    
    
    L<--0.5*no*log(2*pi)-0.5*no*log(sigma2eps)+log_det1-
      (0.5*sigma2eps^{-1})*rotvooinvro
    return(as.numeric(-L))
    
  } # correct efficient negative log likelihood of rho,sigmas- B method
  
  
  # par,x,yo,yotyo,yotBo,BotBo,w,no,I,Bo,model,wpluswt,wtw
  
  
  log_like_full<-function(para,X,y,I,w,no,model){
    # no=ns
    # para<-c(1,2,0.1,0.3,0.4)
    x<-X
    xo<-x[1:no,]
    yo<-y[1:no,]
    
    beta<-para[1:ncol(xo)]
    rho<-para[ncol(xo)+1]
    sigma2eps<-para[ncol(xo)+2]
    sigma2y=para[ncol(xo)+3]
    
    theta<-sigma2y/sigma2eps
    
    
    A<-I-rho*w
    AtA=t(A)%*%A
    v=(I+theta*solve(AtA))
    V=forceSymmetric(v) # make symetric
    V_oo=V[1:no,1:no]
    # det(V_oo)
    if(model==1){
      ro<-yo-xo%*%beta
    }else{
      ro<-yo-(solve(A,x)[(1:no),])%*%beta
    }
    
    
    logdet<-determinant((V_oo),logarithm=T)$modulus[[1]]
    
    L<--0.5*no*log(2*pi)-0.5*no*log(sigma2eps)-0.5*logdet-
      (0.5*sigma2eps^{-1})*t(ro)%*%solve(V_oo,ro)
    return(as.numeric(-L))
    
  } # My marginal log-like full variables
  
  
  
  sigma2eps<-omega_hat
  sigma2y<-theta_hat*omega_hat
  
  estimates<-c(as.vector(beta_hat),as.numeric(rho_hat),
               as.numeric(sigma2eps),as.numeric(sigma2y))
  
  
  
  t2<-Sys.time()
  if(Hessian){
    
    if(h.full.numerical){
      hesmat_full<-numDeriv::hessian(func=log_like_full,
                                     x=estimates,
                                     I=I,X=x,y=y,w=w,no=ns,
                                     model=model)
      # diag(solve(hesmat_full))
    }else{
      # print("B method is used when calculating Hessian.....")
      Hess.rho.sigmas<-numDeriv::hessian(func=log_like_B,beta=as.matrix(beta_hat,ncol=1),
                                         x=estimates[(ncol(x)+1):length(estimates)],
                                         I=I,X=x_till,
                                         yo=yo,yotyo=yotyo,yotBo=yotBo,BotBo=BotBo,wpluswt=wpluswt,
                                         wtw=wtw,w=w,no=no,Bo=Bo,method.args=list(eps=1e-4,
                                                                                  d=0.0001,
                                                                                  zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      
      
      nofbetas<-ncol(x_till_o)
      info<-matrix(rep(0,(ncol(x)+3)^2),ncol=ncol(x)+3)
      
      # ro<-y_o-x_till_o%*%beta_hat
      # Bot.ro<-t(Bo)%*%ro
      # rotvooinvro=t(ro)%*%ro-theta*t(Bot.ro)%*%solve(chol.AtAplus,Bot.ro)
      # omega_hat=rotvooinvro/no
      
      x_till_otvooin_x_till_o=beta_1_inner
      exactinfo.beta<-as(x_till_otvooin_x_till_o*as.numeric((1/sigma2eps)),"matrix")
      info[1:nofbetas,1:nofbetas]<-exactinfo.beta
      info[(nofbetas+1):nrow(info),(nofbetas+1):ncol(info)]<-Hess.rho.sigmas
      hesmat_full<-info
      # solve(hesmat_full)
    }
    
    
  }else{
    hesmat_full<-0
    hess_2<-0
  }
  t3<-Sys.time()
  Hess.time<-t3-t2
  Est.time<-t2-t1
  Times<-c(Est.time,Hess.time)
  names(Times)<-c("Est.time","Hess.time")
  
  estimates.list=list(beta=beta_hat,rho=rho_hat,sigma2eps=sigma2eps,sigma2y=sigma2y)
  
  # return(list(Estimates=l,Hessian=hesmat_full,Times=Times))
  info.beta.rho<-0
  return(list(Estimates.vec=estimates,Hessian=hesmat_full,
              Estimates.list=estimates.list,Times=Times,optim.result=optimresult,info.beta.rho=info.beta.rho))
  
} # MML-P algorithm,

#################### Data simulation #####################

simulateSEM<-function(para,p,weightmat){
  
  
  weightmat<-as(weightmat,"CsparseMatrix")
  sigma2=para[p+3]
  rho=para[p+2]
  n=ncol(weightmat)
  m=0
  std=1
  X<-matrix(rep(c(1),n))
  
  for(i in 1:p){
    X<-cbind(X,rnorm(n,m,std))
  }
  
  b<-para[0:p+1]
  errors_part<-rnorm(n,sd=sqrt(sigma2))
  
  Y<-X%*%b+solve(Diagonal(n)-rho*weightmat,errors_part)
  
  return(list("Independent"=X[,2:(p+1)],"Dependent"=Y))
  
} # Generate SEM, efficient generator

simulateSAM<-function(para,p,weightmat){
  
  
  std_mat<-weightmat
  # g<-graph.lattice(c(sqrtn,sqrtn)) 
  # g<-connect.neighborhood(g,1) #connect all vertices by one edge.
  # wmat<-as_adjacency_matrix(g)
  # row_sums <- rowSums(wmat)
  # factors <- 1 / (row_sums)
  # D <- Diagonal(x=factors)
  # std_mat <- D %*% wmat
  
  
  simSAM<-function(para,n,r_stdwm,disx=c(0,1)){
    
    nofbeta<-length(para)-2
    sima2=para[nofbeta+2]
    rho=para[nofbeta+1]
    n=ncol(r_stdwm)
    
    I<-Diagonal(n)
    
    m<-disx[1]
    std<-disx[2]
    
    x<-matrix(rep(c(1),n))
    
    for(i in 1:(nofbeta-1)){
      x<-cbind(x,rnorm(n,m,std))
    }
    
    b<-para[1:nofbeta]
    
    errors_part<-rnorm(n,sd=sqrt(sima2))
    
    sp_part=(I-rho*r_stdwm)
    y1<-x%*%b+errors_part
    
    y<-solve(sp_part,y1)
    x<-as.matrix(x[,-1])
    
    return(list("Independent"=x,"Dependent"=as.matrix(y)))
    
  } # Generate SAM,
  
  n=ncol(std_mat)
  
  syntheticdata<-simSAM(para=para,n=n,r_stdwm=std_mat)
  
  rest<-c(syntheticdata,"std.weighmat"=std_mat)
  
  return(rest)
  
  
  
} # Generate SAM, efficient generator

splitter<-function(x,y,w,wmat,p_sample){
  
  n=nrow(w)
  ns=round(n*p_sample)
  nu=n-ns
  x=as.matrix(x,nrow=n)
  
  s<-sample(1:n,ns) #sampled
  s<-sort(s)
  u<-setdiff(1:n,s)  # non-sample units
  
  ns<-length(s)
  nu<-n-ns
  
  ys=y[s] #Sold unsold y
  yu=y[u]
  
  xs=as.matrix(x[s,],nrow=ns) #sold unsold x
  xu=as.matrix(x[u,],nrow=nu)
  
  reordered_wmat<-wmat[c(s,u),c(s,u)]
  reordered_wmat_oo<-as(reordered_wmat[1:ns,1:ns], "CsparseMatrix")
  
  row_normalize <- function(mat) {
    mat / rowSums(mat)
  }
  reordered_wmat_oo<-row_normalize(reordered_wmat_oo)
  reordered_wmat_oo[is.na(reordered_wmat_oo)] <- 0
  
  # reordered_wmat_oo<-mat2listw(reordered_wmat_oo,style="W")
  # reordered_wmat_oo=listw2mat(reordered_wmat_oo)
  # reordered_wmat_oo<-as(reordered_wmat_oo,"CsparseMatrix")
  
  
  
  w<-w[c(s,u),c(s,u)] # Create w partitioned matrix
  y<-cbind(c(ys,yu)) # create full (observed+unobserved) x
  x<-xs
  x<-rbind(x,xu) # create full (observed+unobserved) y
  
  l=list("W"=w,"Y"=y,"X"=x,"ns"=ns,"w_oo"=reordered_wmat_oo,"y_o"=ys,"x_o"=xs)
  return(l)
  
} #Make missing data, and partitions of y,W,X, accordingly.


#################### Function to summarize results #####################


createCI_1<-function(variances,estimate,truepara){
  #95% confidence interval
  
  errors<-qnorm(0.975)*sqrt(variances)
  
  CI<-matrix(rep(0,7*length(estimate)),ncol=7)
  # i=5
  for(i in 1:length(estimate)){
    lower<-estimate[i]-errors[i]
    upper<-estimate[i]+errors[i]
    CI[i,]<-c(truepara[i],estimate[i],(estimate[i]-truepara[i])^2,variances[i],lower,upper,ifelse((truepara[i]>=lower&&truepara[i]<= upper), 1, 0))
  }
  
  colnames(CI)<-c("true.val","estimates",
                  "squared.error","var",
                  "L_CI","U_CI","cover")
  
  return(CI)
} #calculate CI



