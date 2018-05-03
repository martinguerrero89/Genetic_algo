#Functions

#Centered pearson distance

gpuDist=function(x){
  x=vclMatrix(t(x))
  mx=colMeans(x)
  x2=colMeans(x^2)
  mx2=mx^2
  sx=sqrt(x2-mx2)
  pDist=((crossprod(x,x)/nrow(x))-(mx %o% mx))/sx%o%sx
  return(as.dist(1-pDist))
}

#Distance functions

mydist=function(c) {gpuDist(c)}                     #c is expression matrix

#To use hierchical clustering 
myclust=function(d) {hclust(d,method="average")}    #d is a distance object

mycluster1 <- function(d, k) {                      #d is a distance object and k is number of clusters
  hc=cutree(myclust(d),k=k)
  return(list(cluster=hc))
}

#To use partition around medioids (PAM), This is the default for the the current aplication

mycluster2<- function(c,k){
  return(list(cluster=pam(c,k,cluster.only=TRUE)))
}

#Humming distance

humdist= function(x,y){sum(abs(x-y))}   #Function to calculate Humming distance between two binary vectors

#Function to calculate the centroids of different groups (classes)
kcentroid=function(data,class){
  L=list()
  c=unique(unlist(class))
  for(i in c){
    if(sum(unlist(class)==i)>1){
      x=rowMeans(data[,unlist(class)==i])
      L[[i]]=x
    }else{
      L[[i]]=data[,unlist(class)==i]
    }}
  L=t(do.call(rbind,L))
  return(L)
}

#Function to calculate the distance to the centroids
corrF= function(data,centroid) apply(data,2,cor,y=centroid,method="spearman")

#Given the distance to the centroids classify the samples
classify= function(data,centroid){
  R=corrF(data,centroid)
  scores<-apply(R,2,which.max)
  return(scores)
}

#harmonic mean tends to be robust with high outlayers (robust for overfitting)
#with high penalty on small values
hmean=function(a){1/mean(1/a)}


#Fitness by RMST (Restricted Mean Survival Time, https://www.bmj.com/content/357/bmj.j2250)
cDist= function(x){ #ad-hoc function, x is the RMST
  d=x[order(x)]     
  l=length(d)-1
  c=c(0,1)
  dif=as.numeric()
  for(i in 1:l){
    dif=c(dif,diff(d[c+i]))
  }
  return(hmean(dif)*l)
}

fitness=function(data=survclass){
  score=tryCatch({ 
    t=survival:::survmean(survfit(Surv(OS,vital_status==vit_stat)~clustclass,data=data),rmean=max(data$OS))[[1]][,"*rmean"] #This function calculates the RMST (comes from package Survival) 
    cDist(t)},error=function(e)return(0)) #If cDist cannot be calculated, the difference is set to 0 (no difference between curves)
  return(score)
}


#crossvalidation function based on: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3105299/
#Data is the expression matrix
#flds is a list with the indexes to partition the data in n different folds
#indv is the solution to test, namely, a binary vector to subset the genes to test
#vit_stat is the codification for an event in the clinical data (e.g. "Dead"), is used in the fitness function

#The function returns two fitness: fit1 (the mean silhouette) and fit2 (the ad-hoc function to estimate the differences between curves)
crossvalidation= function(Data, flds,indv,k,vit_stat=vit_stat){
  Data= Data[indv,] 
  D=mydist(t(Data))
  clustclass=NULL
  for(i in 1:length(flds)){
    trainFold=as.vector(unlist(flds[-i]))
    testFold=as.vector(unlist(flds[i]))
    trainData=Data[,trainFold]
    testData= Data[,testFold]
    sub=subset(D,trainFold)
    hc=mycluster2(sub,k)
    CLASS= classify(testData,kcentroid(trainData,hc$cluster))
    clustclass=c(clustclass,CLASS)
  }
  
  clustclass=clustclass[order(as.vector(unlist(flds)))]
  fit1=mean(silhouette(clustclass,D)[,3])
  survclass= cbind(clustclass,clinical[,c("OS","vital_status")])
  fit2= fitness(survclass)
  
  return(c(fit1,fit2))
}

#Modified Binary lector to read cluster codification (3 bit binary)
BinToDec <- function(x) 
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))+2

DecToBin<- function(x){
  x=x-2
  rev(as.integer(intToBits(x))[1:3])
}

#Minimum number of genes to use in a solution (A constraint for the algorithm)
minGenes= function(x){ sum(x)>=4}

#Multiple point crossover
#a and b is solution 1 and 2 respectively (binary vectors)
#n is the number of cut points
crossover=function(a,b,n){
  if(length(a)!=length(b)) {stop("vectors of unequal length")}
  l=length(a)
  if(n>=(length(a)-1)) {stop("number of cut points bigger than possible sites")}
  points=sample(2:(l-1),n,replace=FALSE)
  to=c(points[order(points)][-n],l)
  from= c(1,to[-length(to)]+1)
  cutpoints=list()
  for(i in 1:n){
    cutpoints[[i]]=seq(from[i],to[i])
  }
  achild=as.numeric()
  bchild=as.numeric()
  for(i in 1:n){
    if(i%%2==0){
      achild=c(achild,a[cutpoints[[i]]])
      bchild=c(bchild,b[cutpoints[[i]]])
    }else{
      achild=c(achild,b[cutpoints[[i]]])
      bchild=c(bchild,a[cutpoints[[i]]])
    }}
  return(list(achild,bchild))
}


