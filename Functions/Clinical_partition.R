#Clinical data perparation and partition

#Check order
if(!identical(colnames(prob_matrix),rownames(clinical))){
  stop("expression matrix and clinical data are not adequately named")
}

#For TCGA data
#vit_stat= "Dead"
#clinical$OS=clinical$OS/30

#For metabric
vit_stat= "Died of Disease"

##
stat= clinical$vital_status== vit_stat
stat= as.numeric(stat)

set.seed(441)
cvcox <- cv.CoxBoost(time=clinical$OS,status=stat,x=t(prob_matrix),K=nCV,type="verweij",maxstepno=500)
flds=cvcox$folds
