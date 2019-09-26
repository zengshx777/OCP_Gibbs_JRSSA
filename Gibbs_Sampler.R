#Gibbs Sampler for Bayes LIV
#S: principal strata

#tau.s: random effect for principal strata
#tau.y.1,tau.y.0: random effect for the outcome
#sigma.s: error in principal strata
#sigma.y: error in outcome model
#beta.s:coefficients for X in principal strata
#beta.y.0,beta.y.1:coefficients for X in outcome model
#alpha.0,alpha.1:threshold for ordinal outcome
#Y.latent.0,Y.latent.1: latent variable for proportional odds
#prior for beta(beta_s,beta_y)
#phi : augmented parameter to obtain t-distribution
#gamma.s.0,gamma.s.1: Coefficients for the principal strata

#Parameter to approximate logistic distribution
#v=7.3
#sigma.y=pi^2*(v-2)/(3*v)

#Input Data
#Y_obs
#Z: Values of Instrumental Variables
#T.ind: Treatment Indicator
#X: Covariates Matrix
#C: number of province
#C.ind:Cluster membership for random Effect
#N:Sample size 
#N.cluster:Number in each cluster
#C.list.in: list of which(C.ind==j)
#C.list.1=N.ind intersect with treat.id
#C.list.0=N.ind intersect with control.id
#K: res.level-1

library(truncnorm)
library(MASS)

#Matrix to restore Posterior Sample
alpha.y.1=alpha.y.0=matrix(0,nrow=mcmc.time,ncol=K)
alpha.y.0[1,]=alpha.y.1[1,]=seq(-100,100,length=K)

#Exception
if(K1!=5)
{
  alpha.y.1[,1:(5-K1)]=-Inf
}
if(K0!=5)
{
  alpha.y.0[,1:(5-K0)]=-Inf
}

Y.latent=matrix(0,nrow=mcmc.time,ncol=N)

beta.y.1=matrix(0,nrow=mcmc.time,ncol=ncov+1)
beta.y.0=matrix(0,nrow=mcmc.time,ncol=ncov+1)
beta.s=matrix(0,nrow=mcmc.time,ncol=ncov+1)
S=matrix(0,nrow=mcmc.time,ncol=N)
S[1,treat.id]=Z[treat.id]-10
S[1,control.id]=Z[control.id]+10
gamma.s.1=gamma.s.0=rep(1,mcmc.time)
phi=matrix(1,nrow=mcmc.time,ncol=N)
sigma.s=rep(1,mcmc.time)
v=7.3
sigma.y=sqrt(pi^2*(v-2)/(3*v))
sigma.y=1.6
tau.y.1=tau.y.0=tau.s=matrix(0,nrow=mcmc.time,ncol=C)
tau.y.0.sigma=tau.y.1.sigma=tau.s.sigma=rep(1,mcmc.time)
#Importance Weight
w=numeric(mcmc.time)

#Prior Parameter
Prec.prior.y=diag(rep(0.01,ncov+1))
Prec.prior.s=diag(rep(0.01,ncov+1))
beta.prior.y=rep(0,ncov+1)
beta.prior.s=rep(0,ncov+1)
#Prior variance of gamma
gamma.sigma.prior=10
#Prior gamma parameter for sigma.s
prior.a=0.5;prior.b=0.5
#Prior gamma parameter for tau.y.0.sigma/tau.y.1.sigma/tau.s.sigma
prior.a.tau=0.2;prior.b.tau=0.2

for (t in 2:mcmc.time)
{
  #Updating Y_latent
  threshold.1=c(-Inf,alpha.y.1[t-1,],Inf)
  threshold.0=c(-Inf,alpha.y.0[t-1,],Inf)
  Y.latent[t,treat.id]<-rtruncnorm(N1,mean=X.aug[treat.id,]%*%beta.y.1[t-1,]+
                           gamma.s.1[t-1]*S[t-1,treat.id]+delta.z*Z[treat.id]+
                           tau.y.1[t-1,C.ind[treat.id]],
                         a=threshold.1[Y.obs[treat.id]],
                         b=threshold.1[Y.obs[treat.id]+1],
                         sd=sigma.y/sqrt(phi[t-1,treat.id]))
  
  Y.latent[t,control.id]<-rtruncnorm(N0,mean=X.aug[control.id,]%*%beta.y.0[t-1,]+
                           gamma.s.0[t-1]*S[t-1,control.id]+
                           tau.y.0[t-1,C.ind[control.id]],
                         a=threshold.0[Y.obs[control.id]],
                         b=threshold.0[Y.obs[control.id]+1],
                         sd=sigma.y/sqrt(phi[t-1,control.id]))
  
  #Updating Principal Strata
  prec.1=1/sigma.s[t-1]^2
  for (i in 1:N)
  {
    if (T.ind[i]==1)
    {
      prec.2=gamma.s.1[t-1]^2/sigma.y^2
      mean.1=X.aug[i,]%*%beta.s[t-1,]+tau.s[t-1,C.ind[i]]
      mean.2=(Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t-1,]-tau.y.1[t-1,C.ind[i]]-delta.z*Z[i])/gamma.s.1[t-1]
      S[t,i]=rtruncnorm(1,mean=mean.1*prec.1+mean.2*prec.2,
                      a=-Inf,b=Z[i],
                 sd=sqrt(1/(prec.1+prec.2)))
    }else{
      prec.2=gamma.s.0[t-1]^2/sigma.y^2
      mean.1=X.aug[i,]%*%beta.s[t-1,]+tau.s[t-1,C.ind[i]]
      mean.2=(Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t-1,]-tau.y.0[t-1,C.ind[i]])/gamma.s.0[t-1]
      S[t,i]=rtruncnorm(1,mean=mean.1*prec.1+mean.2*prec.2,
                 a=Z[i],b=Inf,
                 sd=sqrt(1/(prec.1+prec.2)))
    }
  }
  
  #Updating beta.y.1
  temp_Sigma=solve(Prec.prior.y+1/sigma.y^2*
                     t(X.aug[treat.id,])%*%diag(phi[t-1,treat.id])%*%X.aug[treat.id,])
  temp_mean=temp_Sigma%*%(Prec.prior.y%*%beta.prior.y+
        1/sigma.y^2*t(X.aug[treat.id,])%*%diag(phi[t-1,treat.id])%*%
    (Y.latent[t,treat.id]-gamma.s.1[t-1]*S[t,treat.id]-tau.y.1[t-1,C.ind[treat.id]]-delta.z*Z[treat.id]))
  beta.y.1[t,]=mvrnorm(n=1,mu=temp_mean,Sigma=temp_Sigma)
  
  #Updating beta.y.0
  temp_Sigma=solve(Prec.prior.y+1/sigma.y^2*
                     t(X.aug[control.id,])%*%diag(phi[t-1,control.id])%*%X.aug[control.id,])
  temp_mean=temp_Sigma%*%(Prec.prior.y%*%beta.prior.y+
    1/sigma.y^2*t(X.aug[control.id,])%*%diag(phi[t-1,control.id])%*%
    (Y.latent[t,control.id]-gamma.s.0[t-1]*S[t,control.id]-tau.y.0[t-1,C.ind[control.id]]))
  beta.y.0[t,]=mvrnorm(n=1,mu=temp_mean,Sigma=temp_Sigma)
  
  #Updating beta.s
  temp_Sigma=solve(Prec.prior.s+1/sigma.s[t-1]^2*
                     t(X.aug)%*%X.aug)
  temp_mean=temp_Sigma%*%(Prec.prior.s%*%beta.prior.s+
        1/sigma.s[t-1]^2*t(X.aug)%*%(S[t,]-tau.s[t-1,C.ind]) )
  beta.s[t,]=mvrnorm(n=1,mu=temp_mean,Sigma=temp_Sigma)
  
  #Updating gamma.s.1 
  prec.1=1/gamma.sigma.prior^2
  prec.2=sum(S[t,treat.id]^2*phi[t-1,treat.id])/sigma.y^2
  gamma.mean=1/sigma.y^2*sum((Y.latent[t,treat.id]-X.aug[treat.id,]%*%beta.y.1[t,]-
                    tau.y.1[t-1,C.ind[treat.id]])*phi[t-1,treat.id]*S[t,treat.id])/(prec.1+prec.2)
  gamma.s.1[t]<-rnorm(1,mean=gamma.mean,sd=sqrt(1/(prec.1+prec.2)))
  
  #Updating gamma.s.0
  prec.2=sum(S[t,control.id]^2*phi[t-1,control.id])/sigma.y^2
  gamma.mean=1/sigma.y^2*sum((Y.latent[t,control.id]-X.aug[control.id,]%*%beta.y.0[t,]-
                                tau.y.0[t-1,C.ind[control.id]])*phi[t-1,control.id]*S[t,control.id])/(prec.1+prec.2)
  gamma.s.0[t]<-rnorm(1,mean=gamma.mean,sd=sqrt(1/(prec.1+prec.2)))
              
  #Updating Random Effect
  for (j in 1:C)
  {
    #Update Random Effect in Principal Strata Model
    prec.1=1/tau.s.sigma[t-1]^2
    prec.2=N.cluster[j]/sigma.s^2
    temp.mean=mean(S[t,C.list[[j]]]-X.aug[C.list[[j]],]%*%beta.s[t,])*prec.2/(prec.1+prec.2)
    
    tau.s[t,j]=rnorm(1,mean=temp.mean,sd=sqrt(1/(prec.1+prec.2)))
    #Update Random Effect in Outcome Model
    prec.1=1/tau.y.0.sigma[t-1]^2
    prec.2=length(C.list.0[[j]])/sigma.y^2
    temp.mean=mean((Y.latent[t,C.list.0[[j]]]-gamma.s.0[t]*S[t,C.list.0[[j]]]-
                      X.aug[C.list.0[[j]],]%*%beta.y.0[t,])*phi[t-1,C.list.0[[j]]])*prec.2/(prec.1+prec.2)
    tau.y.0[t,j]=rnorm(1,mean=temp.mean,sd=sqrt(1/(prec.1+prec.2)))
    
    prec.1=1/tau.y.1.sigma[t-1]^2
    prec.2=length(C.list.1[[j]])/sigma.y^2
    temp.mean=mean((Y.latent[t,C.list.1[[j]]]-gamma.s.1[t]*S[t,C.list.1[[j]]]-delta.z*Z[C.list.1[[j]]]-
                      X.aug[C.list.1[[j]],]%*%beta.y.1[t,])*phi[t-1,C.list.1[[j]]])*prec.2/(prec.1+prec.2)
    tau.y.1[t,j]=rnorm(1,mean=temp.mean,sd=sqrt(1/(prec.1+prec.2)))
  }
  
  #Updating Variance for Random Effect
  a=prior.a.tau+C/2
  b=prior.b.tau+0.5*sum(tau.s[t,]^2)
  tau.s.sigma[t]=sqrt(1/rgamma(1,shape=a,rate=b))
  
  b=prior.b.tau+0.5*sum(tau.y.1[t,]^2)
  tau.y.1.sigma[t]=sqrt(1/rgamma(1,shape=a,rate=b))
  
  b=prior.b.tau+0.5*sum(tau.y.0[t,]^2)
  tau.y.0.sigma[t]=sqrt(1/rgamma(1,shape=a,rate=b))
  
  #Updating sigma.s
  a=prior.a+N/2
  b=prior.b+0.5*sum ((S[t,]-X.aug%*%beta.s[t,]-tau.s[t,C.ind])^2)
  sigma.s[t]=sqrt(1/rgamma(1,shape=a,rate=b))
  
  #Updating phi
  if(t_approx==1){
  a=(v+1)/2
  for (i in 1:N)
  {
    if(T.ind[i]==1)
    {
    b=0.5*(v+(Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t,]-
         tau.y.1[t,C.ind[i]]-gamma.s.1[t]*S[t,i]-delta.z*Z[i])^2)/sigma.y^2
    phi[t,i]=rgamma(1,shape=a,rate=b)
    if(phi[t,i]<0.01)
      phi[t,i]=0.01
    }else
    {
      b=v/2+0.5*(Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t,]-
           tau.y.0[t,C.ind[i]]-gamma.s.0[t]*S[t,i])^2/sigma.y^2
      phi[t,i]=rgamma(1,shape=a,rate=b)
      if(phi[t,i]<0.01)
        phi[t,i]=0.01
    }
  }
  }


  
  #Updating Threshold
  #K=4
  for (k in (6-K1):K)
  {
    alpha.y.1[t,k]=runif(1,max(Y.latent[t,intersect(treat.id,which(Y.obs==k))]),
                       min(Y.latent[t,intersect(treat.id,which(Y.obs==k+1))]))
  }
  for (k in (6-K0):K)
  {
    alpha.y.0[t,k]=runif(1,max(Y.latent[t,intersect(control.id,which(Y.obs==k))]),
                       min(Y.latent[t,intersect(control.id,which(Y.obs==k+1))]))
  }
  
  #Calculate Weight
  if(t_approx==1){
  for (i in 1:N)
  {
    if(T.ind[i]==1)
    {
      w[t]=w[t]+dlogis(Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t,]-tau.y.1[t,C.ind[i]]-gamma.s.1[t]*S[t,i]-delta.z*Z[i], 
                 scale = 1, log =T)-
        dt((Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t,]-tau.y.1[t,C.ind[i]]-gamma.s.1[t]*S[t,i])/sigma.y,
            df=v, ncp=0,log =T)
        # dnorm(Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t,]-tau.y.1[t,C.ind[i]]-gamma.s.1[t]*S[t,i]-delta.z*Z[i],sd=sigma.y/sqrt(phi[t,i]),
        #    log =T)
    }else
    {
      w[t]=w[t]+dlogis(Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t,]-tau.y.0[t,C.ind[i]]-gamma.s.0[t]*S[t,i], 
                       scale = 1, log =T)-
        dt((Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t,]-tau.y.0[t,C.ind[i]]-gamma.s.0[t]*S[t,i])/sigma.y,
           df=v, ncp=0,log =T)
      # dnorm(Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t,]-tau.y.0[t,C.ind[i]]-gamma.s.0[t]*S[t,i],sd=sigma.y/sqrt(phi[t,i]),
      #     log =T)
    }
  }}else{
    for (i in 1:N)
    {
      if(T.ind[i]==1)
      {
        w[t]=w[t]+dlogis(Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t,]-tau.y.1[t,C.ind[i]]-gamma.s.1[t]*S[t,i]-delta.z*Z[i], 
                         scale = 1, log =T)-
          # dt((Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t,]-tau.y.1[t,C.ind[i]]-gamma.s.1[t]*S[t,i])/sigma.y,
          #    df=v, ncp=0,log =T)
        dnorm(Y.latent[t,i]-X.aug[i,]%*%beta.y.1[t,]-tau.y.1[t,C.ind[i]]-gamma.s.1[t]*S[t,i]-delta.z*Z[i],sd=sigma.y/sqrt(phi[t,i]),
           log =T)
      }else
      {
        w[t]=w[t]+dlogis(Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t,]-tau.y.0[t,C.ind[i]]-gamma.s.0[t]*S[t,i], 
                         scale = 1, log =T)-
          # dt((Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t,]-tau.y.0[t,C.ind[i]]-gamma.s.0[t]*S[t,i])/sigma.y,
          #    df=v, ncp=0,log =T)
        dnorm(Y.latent[t,i]-X.aug[i,]%*%beta.y.0[t,]-tau.y.0[t,C.ind[i]]-gamma.s.0[t]*S[t,i],sd=sigma.y/sqrt(phi[t,i]),
            log =T)
      }
    }
  }
  if(t%%(mcmc.time/100)==0){
  print(paste("MCMC==",t," th Iterations==",t/mcmc.time*100,"%finished"))
  }
}
mcmc.ind=seq(mcmc.time*6/10,mcmc.time,by=mcmc.time/1000)[-1]
#plot(exp(w[mcmc.ind])/mean(exp(w[mcmc.ind])),type='l')
w_adj=w[mcmc.ind]-mean(w[mcmc.ind])
w_adjust=exp(w_adj)/mean(exp(w_adj))
#sum(0.01<w_adjust&&w_adjust<100)

#Use the empirical distribution of X and Threshold Value
#Extract Threshold Values
y1_mean=beta.y.1[mcmc.ind,]%*%t(X.aug)+tau.y.1[mcmc.ind,C.ind]
y1_mean=t(t(y1_mean)+delta.z*Z)
y0_mean=beta.y.0[mcmc.ind,]%*%t(X.aug)+tau.y.0[mcmc.ind,C.ind]
#Obtain the individual with principal strata within the range of IV
#Complier Index
c.id<-apply(S[mcmc.ind,],1,FUN=function(x){which(zmin<x&x<zmax)})
Y_1_latent<-lapply(1:length(mcmc.ind),FUN=function(x){y1_mean[x,c.id[[x]]]+gamma.s.1[mcmc.ind[x]]*S[mcmc.ind[x],c.id[[x]]]})
Y_0_latent<-lapply(1:length(mcmc.ind),FUN=function(x){y0_mean[x,c.id[[x]]]+gamma.s.0[mcmc.ind[x]]*S[mcmc.ind[x],c.id[[x]]]})

logit<-function(x){1/(1+exp(-x))}
Y_1_expect=lapply(1:length(mcmc.ind),FUN=function(x)
{
  unlist(lapply(1:length(c.id[[x]]),FUN=function(y)
  {
    #K+1-sum(logit(alpha.y.1[mcmc.ind[x],]-y))
    K+1-sum(pnorm((alpha.y.1[mcmc.ind[x],]-Y_1_latent[[x]][y])*sqrt(phi[mcmc.ind[x],c.id[[x]][y]])/sigma.y))
  }))
})

Y_0_expect=lapply(1:length(mcmc.ind),FUN=function(x)
{
  unlist(lapply(1:length(c.id[[x]]),FUN=function(y)
  {
  #  K+1-sum(logit(alpha.y.0[mcmc.ind[x],]-y))
    K+1-sum(pnorm((alpha.y.0[mcmc.ind[x],]-Y_0_latent[[x]][y])*sqrt(phi[mcmc.ind[x],c.id[[x]][y]])/sigma.y))
  }))
})
prte=unlist(lapply(Y_1_expect,mean))-unlist(lapply(Y_0_expect,mean))
###Calculating PRTE
mean(prte)
quantile(prte,0.025)
quantile(prte,0.975)

###Calculating ATT
Y_1_latent<-lapply(1:length(mcmc.ind),FUN=function(x){y1_mean[x,treat.id]+gamma.s.1[mcmc.ind[x]]*S[mcmc.ind[x],treat.id]})
Y_0_latent<-lapply(1:length(mcmc.ind),FUN=function(x){y0_mean[x,treat.id]+gamma.s.0[mcmc.ind[x]]*S[mcmc.ind[x],treat.id]})

Y_1_expect=lapply(1:length(mcmc.ind),FUN=function(x)
{
  unlist(lapply(1:length(treat.id),FUN=function(y)
  {
    #K+1-sum(logit(alpha.y.1[mcmc.ind[x],]-y))
    K+1-sum(pnorm((alpha.y.1[mcmc.ind[x],]-Y_1_latent[[x]][y])*sqrt(phi[mcmc.ind[x],treat.id[y]])/sigma.y))
  }))
})

Y_0_expect=lapply(1:length(mcmc.ind),FUN=function(x)
{
  unlist(lapply(1:length(treat.id),FUN=function(y)
  {
    #  K+1-sum(logit(alpha.y.0[mcmc.ind[x],]-y))
    K+1-sum(pnorm((alpha.y.0[mcmc.ind[x],]-Y_0_latent[[x]][y])*sqrt(phi[mcmc.ind[x],treat.id[y]])/sigma.y))
  }))
})

# logit<-function(x){1/(1+exp(-x))}
# Y_1_expect=lapply(1:length(mcmc.ind),FUN=function(x)
# {
#   unlist(lapply(Y_1_latent[[x]],FUN=function(y)
#   {
#     K+1-sum(logit(alpha.y.1[mcmc.ind[x],]-y))
#   }))
# })
# 
# Y_0_expect=lapply(1:length(mcmc.ind),FUN=function(x)
# {
#   unlist(lapply(Y_0_latent[[x]],FUN=function(y)
#   {
#     K+1-sum(logit(alpha.y.0[mcmc.ind[x],]-y))
#   }))
# })
att=unlist(lapply(Y_1_expect,mean))-unlist(lapply(Y_0_expect,mean))
mean(att)
quantile(att,0.025)
quantile(att,0.975)
