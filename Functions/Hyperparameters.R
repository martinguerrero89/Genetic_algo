#Hyperparameters

population= 50                   #Number of individuals to evaluate
elite=5                          #Number of individuals to pass to next generation unaltered
generations=4                   #Number of generations
chrom_length= length(HSPvector)   #length of chromosome
nCV=5                             #Number of crossvalidations for function "crossvalidation"
COpoints= round(chrom_length/10)  #number of crossing over points
plotgr=FALSE                      #set to TRUE if you want a generation plot
GPU= TRUE                         # to use gpuR


