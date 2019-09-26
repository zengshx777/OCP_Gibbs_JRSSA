#Preprocessing
temp.data=data.jags[[sub_group]]
X=as.matrix(temp.data[,covariate_index])
X.aug=cbind(1,X)
Z=iend-temp.data$ifppr
zmin=min(Z)
zmax=max(Z)
T.ind=temp.data$onechild
K=length(unique(temp.data[,response_index[res.id]]))-1
Y.obs=temp.data[,response_index[res.id]]
ncov=length(covariate_index)

N=nrow(temp.data)
C=length(table(temp.data$province))
N.cluster=table(temp.data$province)
C.ind=temp.data$rd.ind
C.list=list(length=C)
C.list.0=list(length=C)
C.list.1=list(length=C)

treat.id=which(temp.data$onechild==1)
control.id=which(temp.data$onechild==0)
K1=length(unique(temp.data[treat.id,response_index[res.id]]))
K0=length(unique(temp.data[treat.id,response_index[res.id]]))
N1=length(treat.id);N0=length(control.id)


for (c in 1:C)
{
  C.list[[c]]=which(C.ind==c)
  C.list.0[[c]]=intersect(which(C.ind==c),control.id)
  C.list.1[[c]]=intersect(which(C.ind==c),treat.id)
}


