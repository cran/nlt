fwtnpperm<-function(input,f,nkeep=2,intercept=TRUE,initboundhandl="reflect",updateboundhandl="add",neighbours=1,closest=FALSE,LocalPred=AdaptPred,mod=sample(1:length(input),(length(input)-nkeep),FALSE)){

#does the 1D single coefficient lifting transform
#inputtype determines which procedure to follow in order to apply the transform
#input is either a vector of points, or a vector of interval lengths, together
#with a startpoint in its first entry. If startpoint is NA, default is used
#(input matches inputtype).
#nkeep is no. points to keep after transform is performed (at least 2).
#intercept is boolean, denoting whether regression with  intercept is required
#boundaryhandling determines which boundary policy to use.
#closest indicates the way to choose neighbours (boolean)
#neighbours is the number of neighbours on each side, or total number if 
#closest is chosen.
#mod is a vector which gives the order of point removal, i.e. remove 5th,3rd,1st,.. 
#observations in this order; it has length (n-nkeep)


X<-input;
I<-intervals(X,initboundhandl)$i;  #creates interval endpoints with
				 	   #edge correction

lengths<-lengthintervals(X,I,type="midpoints",neighbours,closest)$lengths;
				    #^^^^^^^^^^^^^^^^^ 
origlengths<-lengths #they are in order(x)

X<-as.row(X);
f<-as.row(f);
nkeep<-max(nkeep,2);  #ensures not too many steps are performed

n<-length(X);

removelist<-NULL;	#contains removed points
lengthsremove<-NULL;    #contains interval lengths of removed points
neighbrs<-list();       #list containing the neighbours of the removed 
			#point at each step

gamlist<-list()
predlist<-NULL
alphalist<-list()
lengthlist<-list()
schemehist<-NULL       	#records the history of prediction method (linear,..)
interhist<-NULL    	#records the history of intercept (T,F)
detailslist<-list()
clolist<-NULL
minindexlist<-NULL
mindetailslist<-list()
indiceslist<-list()
history<-list()

pointsin<-matrix(1:n,1,n);
pointsin<-pointsin[order(X)];

coeff<-f;

for (j in 1:(n-nkeep)){

remove<-mod[j];   #removes from list at jth step the element given by 
		  #the jth position of vec (if mod[1]=5, then we remove the 5th observed point)

removelist[j]<-remove;       #updates list of removed points

out<-getnbrs(X,remove,pointsin,neighbours,closest);

nbrs<-out$n;

index<-out$index;

res<-LocalPred(pointsin,X,coeff,nbrs,remove,intercept,neighbours)

if(length(res)==2){
l<-res[[1]]
newinfo<-res[[2]]
nbrs<-newinfo[[3]]
index<-newinfo[[4]]
clolist[j]<-newinfo[[1]]
minindexlist[j]<-newinfo[[2]]
mindetailslist[[j]]<-newinfo[[5]]
indiceslist[[j]]<-newinfo[[6]]
history[[j]]<-c(indiceslist[[j]][[minindexlist[j]]],minindexlist[j])
		}
else{l<-res}

neighbrs[[j]]<-nbrs  #puts appropriate nbrs into neighbrs according to 
		     #LocalPred choice
weights<-l[[4]];
pred<-l[[5]];
bhat<-l[[3]];

if (length(l) == 6) {
            scheme <- NULL
            int <- NULL
            details <- NULL
        }
else{
scheme<-l[[8]];
int<-l[[7]];
details<-l[[9]]
}
coeff[remove]<-coeff[remove]-pred

l1<-PointsUpdate(X,coeff,nbrs,index,remove,pointsin,weights,lengths,updateboundhandl);
coeff<-l1$coeff;
lengths<-l1$lengths;
r<-l1$r;
weights<-l1$weights;
N<-l1$N;
alpha<-l1$alpha;
		
lengthsremove[j]<-lengths[r];
gamlist[[j]]<-weights
predlist[j]<-pred
alphalist[[j]]<-alpha
schemehist[j]<-scheme
interhist[j]<-int
detailslist[[j]]<-details

lengths<-lengths[setdiff(1:length(pointsin),r)];
pointsin<-setdiff(pointsin,remove);

lengthlist[[j]]<-lengths
}  #end j for loop

N<-length(pointsin);

return(list(X=X,coeff=coeff,origlengths=origlengths,lengths=lengths,lengthsremove=lengthsremove,pointsin=pointsin,removelist=removelist,neighbrs=neighbrs,neighbours=neighbours,gamlist=gamlist,alphalist=alphalist,schemehist=schemehist,interhist=interhist,clolist=clolist));
}







