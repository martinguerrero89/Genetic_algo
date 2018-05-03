source("./Functions/Functions.R")
source("./Functions/libraries.R")

clinicalV= read.table("./Data/clinicalValidation.txt",sep="\t")
prob_matrixV= read.table("./Data/validationSet.txt",sep="\t")

resultdir="./Results/"

load(paste0(resultdir,"results.rda"))
C=read.table(paste0(resultdir,"centroids.txt"))
R=rownames(C)
k=PVAL$k[max(PVAL$Generation)]


#Calcular centroides del dataset  completo

#hsp_class=mycluster(t(prob_matrixV[R,]),k)

#C=

#write.table(C,paste0(resultdir,"centroids.txt"),sep="\t",col.names =LETTERS[1:k],row.names=TRUE)
hsp_class= classify(prob_matrixV[R,],C)
##

hsp_classdf= as.data.frame(hsp_class)
hspClinic= cbind(clinicalV,hsp_classdf)

#
#colnames(hspClinic)[21]="clustclass"
#colnames(hspClinic)[21]="hsp_class"
#
vit_stat="Died of Disease"

mysurv <- Surv(hspClinic$OS, hspClinic$vital_status==vit_stat)
tumortotal <- survfit(mysurv~ hspClinic$hsp_class)
totalsdf <- survdiff(mysurv~ hspClinic$hsp_class)
tumortotalpval <- 1 - pchisq(totalsdf$chisq, length(totalsdf$n) - 1)
tumortotalpval <- format(tumortotalpval, digits=4)
COLS= levels(as.factor(hspClinic[,"hsp_class"]))
#pdf(paste0(resultdir,"GA_Kaplan.pdf"))
par(cex=1.35, mar=c(3.8, 3.8, 2.5, 2.5) + 0.1)
plot(tumortotal, main="Gene_class Survival", yscale=100,conf.int=FALSE, col=COLS, lty=1, lwd=2, cex=0.8, xlab="time (months)", ylab="survival(%)", xlim=c(0,max(hspClinic$OS)+0.1*max(hspClinic$OS)));legend("topright", legend=COLS,lty=c(1:1), col=COLS,cex=0.95, bty="n");legend("bottomleft", legend = c("p =", tumortotalpval), horiz=TRUE, xjust=0, bty="n", cex=0.95, x.intersp=0)
#dev.off()

pdf(paste0(resultdir,"expression_heatmap.pdf"))
heatmap.2(as.matrix(prob_matrix[R,]),hclustfun=myclust, distfun=mydist,trace="none",col=redgreen,breaks=seq(-3,3,by=0.05))
dev.off()


pdf(paste0(resultdir,"centroids.pdf"))
heatmap.2(C,col=redgreen,breaks=seq(-0.5,0.5,by=0.05))
dev.off()

#Repetido?
#pdf(paste0(resultdir,"RESULT.pdf"))
#par(cex=1.35, mar=c(3.8, 3.8, 2.5, 2.5) + 0.1)
#plot(tumortotal, main="HSP_class Survival", yscale=100, conf.int=FALSE, col=COLS, lty=1, lwd=2, cex=0.8, xlab="time (months)", ylab="survival(%)", xlim=c(0,max(hspClinic$OS)+0.1*max(hspClinic$OS)));legend("topright", legend=COLS,lty=c(1:1), col=COLS,cex=0.95, bty="n");legend("bottomleft", legend = c("p =", tumortotalpval), horiz=TRUE, xjust=0, bty="n", cex=0.95, x.intersp=0)
#dev.off()

pdf(paste0(resultdir,"Fitness.pdf"))
plot(PVAL[,2],type="l",col="red",lwd=4,xlab="Generation",ylab="Fitness",main="Fitness evolution")
dev.off()

freq=colSums(do.call(rbind,PVAL[,"Genes"]))
index=order(freq,decreasing=T)
count=freq[index]
count[count==0]=NA
#gene names
df <- data.frame(HSP= HSPvector[index], count = count)
#gene index
#df <- data.frame(HSP= c(1:length(HSPvector))[index], count = freq[index])

df$HSP <- factor(df$HSP, levels=df$HSP)
df=na.omit(df)
g=ggplot(df, aes(x=HSP, y=count,fill=-1*count)) + geom_bar(stat="identity",width=1)
g=g+theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))
pdf(paste0(resultdir,"wordfreq.pdf",height=3.5,width=10))
print(g)
dev.off()
#Word cloud

set.seed(1234)
pdf(paste0(resultdir,"wordcloud.pdf"))
wordcloud(words = as.character(df$HSP), freq = df$count, min.freq = 1,
          max.words=100, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(11, "Spectral"))
dev.off()


#to try cv homogeneity
class=as.numeric()
for(i in 1:length(flds)){
  class[flds[[i]]]=i
}
hsp_classcv=class

mysurv <- Surv(hspClinic$OS, hspClinic$vital_status==vit_stat)
tumortotal <- survfit(mysurv~ hsp_classcv)
totalsdf <- survdiff(mysurv~ hsp_classcv)
tumortotalpval <- 1 - pchisq(totalsdf$chisq, length(totalsdf$n) - 1)
tumortotalpval <- format(tumortotalpval, digits=4)
COLS= levels(as.factor(hsp_classcv))
pdf(paste0(resultdir,"CV_kaplan.pdf"))
par(cex=1.35, mar=c(3.8, 3.8, 2.5, 2.5) + 0.1)
plot(tumortotal, main="CV_class Survival", yscale=100, conf.int=FALSE, col=COLS, lty=1, lwd=2, cex=0.8, xlab="time (months)", ylab="survival(%)", xlim=c(0,max(hspClinic$OS)+0.1*(max(hspClinic$OS))));legend("topright", legend=COLS,lty=c(1:1), col=COLS,cex=0.95, bty="n");legend("bottomleft", legend = c("p =", tumortotalpval), horiz=TRUE, xjust=0, bty="n", cex=0.95, x.intersp=0)
dev.off()

