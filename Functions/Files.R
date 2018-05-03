#Files

clinical= read.table("./Data/clinicalGA.txt",sep="\t")
prob_matrix= read.table("./Data/prob_matrixGA.txt",sep="\t")
HSPvector= read.table("./Data/Gene_Vector.txt",sep="\t")
HSPvector= as.character(unlist(HSPvector))
