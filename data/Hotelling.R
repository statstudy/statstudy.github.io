"Hotelling" <- function(da,mu){
# The data matrix is "da".
# The mean vector is mu. (a list of numbers).
#
if(!is.matrix(da))da=as.matrix(da)
Hotelling=NULL
nr = dim(da)[1]
nc = dim(da)[2]
cm=matrix(colMeans(da),nc,1)
S = cov(da)
si=solve(S)
mu0=matrix(mu,nc,1)
dev=cm-mu0 
T2 = nr*(t(dev)%*%si%*%dev)
d2=nr-nc
tt=T2*d2/(nc*(nr-1))
pvalue=1-pf(tt,nc,d2)
Hotelling=cbind(Hotelling,c(T2,pvalue,nc,d2))
row.names(Hotelling)=c("Hoteliing-T2","p.value","1df","2df")
Hotelling
}

