nlt<-function(x,f,J,Pred=AdaptPred,neighbours=1,closest=FALSE,intercept=TRUE,nkeep=2,trule="median"){

n<-length(x)
vec<-matrix(0,J,n-2)
denf<-list()
ghatmat<-matrix(0,J,n)
ghatnat<-matrix(0,1,n)

for (i in 1:J){
v<-sample(1:n,(n-2),FALSE)
vec[i,]<-as.row(v)
denf[[i]]<-denoiseperm(x,f,pred=Pred,neigh=neighbours,int=intercept,clo=closest,keep=nkeep,rule=trule,per=v)
ghatmat[i,]<-as.row(denf[[i]]$fhat$coeff)
}

aveghat<-apply(ghatmat,2,mean)

df<-denoise(x,f,pred=Pred,neigh=neighbours,int=intercept,clo=closest,keep=nkeep,rule=trule)
ghatnat<-as.row(df$fhat$coeff)



return(list(vec=vec,ghatmat=ghatmat,ghatnat=ghatnat,aveghat=aveghat))

}


