#raffonei-waglerianus

# set working directory
setwd("/Users/mgabrielli/Documents/Postdoc/Projects/Brazil/Belem2025/workshop/")

# read files produced by plink
# eigenvalues
pca.plink.all.eig=read.table("Podarcis.raffonei-waglerianus.allchr.miss0.75.thinned.eigenval")
# eigenvectors
pca.plink.all=read.table("Podarcis.raffonei-waglerianus.allchr.miss0.75.thinned.eigenvec")

# extract individual names and create a matrix with it
ind=pca.plink.all[,1]
ind=as.matrix(ind)
colnames(ind)="ind"

# read metadata information
library(data.table)
mydata<-fread("Podarcis_inds_sp_pop.txt",header=FALSE)
colnames(mydata)=c("ind","sp","pop")

# merge metadata information with individual name, keeping the order of individuals in plink files
tab=merge(ind,mydata,by="ind",sort=FALSE)
pop=factor(tab$pop)

# remove the first 2 columns that contain individual names
pca.plink.all=pca.plink.all[,-c(1,2)]

# plot the first 2 axes of the PCA and write it in a PDF file
pdf("PCA_raffonei_waglerianus.pdf",5,5)
plot(pca.plink.all[,2]~pca.plink.all[,1],pch=16,border="black",col=c("goldenrod","chocolate4","#009999")[as.factor(pop)],ylab="",xaxt='n',yaxt="n",xlab="",main=("Podarcis raffonei / waglerianus"))
legend("topleft",c(expression(paste(italic("P. raffonei")," LC")), expression(paste(italic("P. raffonei")," ST")),expression(paste(italic("P. waglerianus")))),pch=16,col=c("goldenrod","chocolate4","#009999"))
title(xlab=paste0("PC1 (",round(pca.plink.all.eig[1,]/sum(pca.plink.all.eig)*100,2),"%)"),line=0)
title(ylab=paste0("PC2 (",round(pca.plink.all.eig[2,]/sum(pca.plink.all.eig)*100,2),"%)"),line=0)
dev.off()
