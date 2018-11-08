library(copula)
rm(list = ls())
library(huge)
library(earth)
library(MASS)
m<-read.table("clipboard",head=T)

X<-m

size<-dim(X)
p<-size[2]
p
tpre<-matrix(c(1),p,p)
myMat <-matrix(0, p,p)

colnames(X) <- paste("X",1:p,sep="")
rownames(myMat) <- colnames(myMat) <- paste("X",1:p,sep="")

t<-colnames(myMat)

for (i2 in 1:p){
  
  y<-X[,i2]
  x<-X[,-i2]
  
  fit1<-earth(y~.,data=x,degree=2,trace=2)
  a <- fit1$dirs
  b<-	fit1$coeff
  a
  b
  
  b1<-rownames(b)
  b2<-dim(b)
  p2<-b2[1]
  
  
  if(p2==1) {
    myMat==myMat
  }else{
    for(i3 in 2:p2){
      
      k<-a[b1[i3],]
      myMat[t[i2], names(which(k!=0))] <- 1
      
      if( length(names(which(k!=0)))==2){
        a99<-names(which(k!=0))
        myMat[a99[1],a99[2]] <- 1
      }
      
    }#end of i3
  }#end of if 
  
}#i2

prec.mat<-myMat

diag(prec.mat)<-1
isSymmetric(prec.mat)
sym.prec.mat<-matrix(0,p,p)
if(isSymmetric(prec.mat)==TRUE)
{sym.prec.mat<-prec.mat}
if(isSymmetric(prec.mat)==FALSE)
{ for(i4 in 1:p){
  for(i5 in 1:p){
    if(prec.mat[i4,i5]==1 || prec.mat[i5,i4]==1)
    { sym.prec.mat[i4,i5]<-1
    sym.prec.mat[i5,i4]<-1}
  }#end of i5
}#end of i4
}
isSymmetric(sym.prec.mat)

epre<-sym.prec.mat
tpre<-as.matrix(tpre)



TP<-0   
FN<-0
FP<-0
TN<-0

for(i8 in 1:p)
{ 
  for(j8 in 1:p)
  { 
    if(tpre[i8,j8]==epre[i8,j8])
    { 
      if(tpre[i8,j8]==1)
        TP=TP+1 
      else
        TN=TN+1
    }
    else
    { 
      if(tpre[i8,j8]==1)
        FN=FN+1 
      else
        FP=FP+1
    }
  }#end of j8.
}#end of i8.



fdr<-round(FP/(FP+TP),3)

fpr<-round(FP/(FP+TN),3)

pre<-round(TP/(TP+FP),3)

rec<-round(TP/(TP+FN),3)

spe<-round(TN/(TN+FP),3)

acc<-round((TP+TN)/(TP+FP+FN+TN),3)

F<-round(2*(pre*rec)/(pre+rec),3)

MCC<-round((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),3)

write(c(i1,TP,FN,FP,TN,spe,fdr,fpr,pre,rec,acc,F),"interaction-effect-VB.txt",ncol=12,append=T,sep = "\t")







