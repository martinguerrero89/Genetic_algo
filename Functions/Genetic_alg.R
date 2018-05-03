#Create folder for results...
if(plotgr==TRUE){
  if(!file.exists("Results")){
    dir.create("./Results/generations",recursive=TRUE)
  }else{
    if(!file.exists("./Results/generations")){
      dir.create("./Results/generations",recursive=TRUE)
    }
  }
}

#parallele computing
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

#Option to run it in GPU mode or CPU mode
if(GPU==TRUE){
  gpupackage="gpuR"
  computingtype= "using GPU computing"
}else{
  gpupackage="amap"
  mydist<<-function(c) {Dist(c,method="pearson")}
  computingtype="using CPU computing"
}

PARETO=list() #Empty list to save the solutions

#1. Create random population of solutions

#Creating random clusters from 2-9 encoded in 3 bits binary
Nclust=replicate(population,sample(c(1,0),size=3,replace=T))

#Matrix with random TRUE false with uniform distribution, representing solutions to test
X=matrix(NA,nrow=population,ncol=chrom_length)
for(i in 1:population){
  prob=runif(1,0,1)
  X[i,]= sample(c(1,0),chrom_length,replace=T,prob=c(prob,1-prob))
}


#####Main loop
for(g in 1:generations){

  start_time <- Sys.time() #Measures generation time
  
  #2.Calculate the fitness f(x) of each chromosome x in the population.
  Fit1= apply(X,1,minGenes) #Apply constraints (min 4 genes per solution)
  X= X[Fit1,]               
  X=apply(X,2,as.logical)
  n=nrow(X)
  Nclust=Nclust[,Fit1]
  
  k= apply(Nclust,2,BinToDec) #Decode binary vector to get number of clusters
  
  
  Fit2=foreach(i=1:nrow(X),.packages=c('cluster',gpupackage,"cba","survival"),.combine=rbind) %dopar% {
    crossvalidation(prob_matrix,flds,X[i,],k[i],vit_stat=vit_stat)
  }   #Calculate Fitnes 1 (silhouette) and 2 (Survival differences)
  
  
  PARETO[[g]]= Fit2 #Saves the fitnes of the solutions of the current generation 
  
  ranking=fastNonDominatedSorting(Fit2*-1) #NonDominatedSorting from package nsga2R. -1 multiplication to alter the sign of the fitness (function sorts from min to max)
  rnkIndex <- integer(n)
  i <- 1
  while (i <= length(ranking)) { #saves the rank of each solution
    rnkIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  
  X1=cbind(X,k,Fit2,rnkIndex) #data.frame with solution vector, number of clusters and ranking
  objRange=apply(Fit2,2,max)-apply(Fit2,2,min) # Range of fitness of the solutions 
  CrowD= crowdingDist4frnt(X1,ranking,objRange) #Crowding distance of each front (nsga2R package)
  CrowD=apply(CrowD,1,sum)
  
  X1= cbind(X1,CrowD) #data.frame with solution vector, number of clusters, ranking and crowding distance
  
  #output for the generation
  print(paste0("Generation ",g," Non-dominated solutions:"))  
  print(X1[X1[,"rnkIndex"]==1,(chrom_length+1):(chrom_length+5)])
  
  
  #Termination parameter
  if(g==generations) stop("the algorithm has ended")
  
  
  #3.Repeat the following steps until n offspring have been created:
  
  Fitesst=order(X1[,"rnkIndex"],X1[,"CrowD"]*-1) #Order solutions according to their ranking and their crowding distance
  
  Y=X1[Fitesst,1:chrom_length]
  Yk= Nclust[,Fitesst]
  
  New= matrix(NA,ncol=chrom_length,nrow=population) #Create empty matrix to add new individuals
  New[1:elite,]= Y[1:elite,]                        #Save elite solutions
  NewK= matrix(NA,nrow=3,ncol=population)           #same for cluster chromosome
  NewK[,1:elite]=Yk[,1:elite]
  
  
  matingPool <- tournamentSelection(X1,population-(elite),6) #Use tournament selection (Size 6), to select perents that will give offspring
  
  
  count=(sum(!is.na(rowSums(New)))) #Count how many offsprings are still needed to reach the original population size
  while(sum(is.na(rowSums(New))) > 0){
    count=count+1
    ##a.Select a pair of parent chromosomes from the matingPool
    Pair=sample(1: nrow(matingPool),2,replace=F)
    
    
    ##b.With probability pc (the "crossover probability" or "crossover rate"), cross over the pair at a
    #n randomly chosen points (chosen with uniform probability) to form two offspring. If no
    #crossover takes place, form two offspring that are exact copies of their respective parents.
    Cp=0.7
    if(sample(c(1,0),1,p=c(Cp,1-Cp))== 1){
      #multiple point crossingover    
      offspring=crossover(matingPool[Pair[1],1:chrom_length],matingPool[Pair[2],1:chrom_length],COpoints)
      off1= offspring[[1]]
      off2= offspring[[2]]
    }else{
      off1= matingPool[Pair[1],1:chrom_length]
      off2= matingPool[Pair[2],1:chrom_length]
    }
    
    #Crossover for cluster N
    clustpool= sapply(matingPool[Pair,"k"],DecToBin)
    Cp=0.5
    if(sample(c(1,0),1,p=c(Cp,1-Cp))== 1){
      point= sample(1:3,1)
      x1= clustpool[1:point-1,1]
      x2= clustpool[point:3,1]
      y1= clustpool[1:point-1,2]
      y2= clustpool[point:3,2]
      offk1= c(x1,y2)
      offk2= c(y1,x2)
    }else{
      offk1= clustpool[,1]
      offk2= clustpool[,2]
    }
    
    
    ##c.Mutate the two offspring at each locus with probability Mp (the mutation probability or mutation rate),
    # and place the resulting chromosomes in the new population.
    #in this case the mutation probability is proportionally inverse to the humming distance
    HD= humdist(matingPool[Pair[1],1:chrom_length],matingPool[Pair[2],1:chrom_length])
    Mp= (chrom_length-HD)/chrom_length
    
    
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      mutpoint= sample(1:chrom_length,1)
      off1[mutpoint]=abs(off1[mutpoint]-1)
    }
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      mutpoint= sample(1:chrom_length,1)
      off2[mutpoint]=abs(off2[mutpoint]-1)
    }
    
    
    #Mutation for cluster N
    Mp=0.01
    
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      kmutpoint= sample(1:3,1)
      offk1[kmutpoint]=abs(offk1[kmutpoint]-1)
    }
    if(sample(c(1,0),1,p=c(Mp,1-Mp))== 1){
      kmutpoint= sample(1:3,1)
      offk2[kmutpoint]=abs(offk2[kmutpoint]-1)
    }
    
    
    #Add offsprings to new generation
    New[count,]= off1
    New[count+1,]= off2
    NewK[,count]=offk1
    NewK[,count+1]=offk2
    
    
  }
  #4.Replace the current population with the new population.
  X=New
  Nclust=NewK
  #5.Go to step 2
  end_time <- Sys.time()
  t=end_time-start_time
  
  print(t)
  gc()
  
  
  
  if(plotgr==TRUE){
    if(g==1){PF=PARETO[[g]]}else{PF=PARETO[[g-1]]}
    png(paste0("./Results/generations/",g,"Pareto.png"), pointsize = 8,width = 840*(population/100), height = 840)
    plot(PF,xlim=c(0,0.6),ylim=c(0,60),cex=3)
    points(matingPool[,(chrom_length+2):(chrom_length+3)],col="red",cex=3,pch=13)
    points(X1[Fitesst[1:elite],(chrom_length+2):(chrom_length+3)],col="green",cex=3,pch=9)
    dev.off()
  }
  
}

stopCluster(cluster)

