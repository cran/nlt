"transmatdualperm" <-
function(x,f,Pred=AdaptPred,neigh=1,int=TRUE,clo=FALSE,keep=2,perm=sample(1:length(x),(length(x)-keep),FALSE)){

#if x has multiple points, then uses fwtnpmp
if (length(x)==length(unique(x))){
out<-fwtnpperm(x,f,LocalPred=Pred,neighbours=neigh,intercept=TRUE,closest=clo,nkeep=keep,mod=perm)
}
else{
out<-fwtnpmp(x,f,LocalPred=Pred,neighbours=neigh,intercept=TRUE,closest=clo,nkeep=keep,mpdet="min")
}

n<-length(out$X)

Alistdual<-list()
Tlistdual<-list()
matno<-n-keep

pointsin<-out$pointsin
remlist<-out$removelist
newpoints<-(c(pointsin,rev(remlist)))

for (i in 1:matno){
Ai<-Amatdual(i,pointsin,remlist,out$neighbrs[[i]],out$gamlist[[i]],out$alphalist[[i]])
Alistdual[[i]]<-Ai$A
}

Tlistdual[[1]]<-Alistdual[[matno]]

for (j in 2:matno){
augment<-rbind(cbind(Tlistdual[[j-1]],0),0)
augment[nrow(augment),nrow(augment)]<-1

Tlistdual[[j]]<-augment%*%Alistdual[[matno-j+1]]

}  #end j


W<-Tlistdual[[matno]]
x<-out$X

Wnew<-matrix(0,length(x),length(x))
reo<-NULL
for (i in 1:length(x)){
reo[i]<-which(x[i]==x[newpoints])
}
for (i in 1:length(x)){
Wnew[,i]<-W[,reo[i]][reo]
		}

#########!!!!Wnew is the matrix of T IN THE ORDER of COEFF, i.e. each column corr to f_1 
#..f_n (in the order of being observed).


return(list(out=out,Wnew=Wnew,x=x))
#return(Alistdual)

}







