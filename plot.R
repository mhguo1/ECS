#NOTE, this is probably the least elegant code you'll ever see, but hopefully it makes sense
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(scales)
require(lattice)

#Generate colors
cols<-brewer.pal(7, "Dark2")
pops<-c("AFR", "AMR", "ASJ", "EAS", "NFE", "SAS", "USA")
cols<-c(cols[1:4], cols[6:7], "gray40")
cols2<-hue_pal()(4)

#read in data
dat<-read.delim("clinvar_subset_summed_table.txt", header=T, stringsAsFactors = F, sep="\t")

#Remove Finnish
dat<-subset(dat, select=-c(FIN, ALL))
names(dat)[1]<-"GENE"

#Calculate max AF across populations
dat[, "MAX_AF"] <- apply(dat[, 2:8], 1, max)

#Subset data by carrier rate at MAX_AF > 0.01, 0.005 and 0.001
dat_01<-subset(dat, MAX_AF>0.01)
dat_001<-subset(dat, MAX_AF>0.001)
dat_005<-subset(dat, MAX_AF>0.005)


#Plot top 10 GCRs by population
#Figure 2A
gcr.all.dat<-data.frame(GENE=character(), GCR=numeric(), GENE=character(), POP_GENE=character(), POP_COLOR=character(), stringsAsFactors = F)
for(p in 1:length(pops)){
  temp<-dat[,c("GENE",pops[p])]
  temp$POP_GENE<-paste(pops[p], temp$GENE, sep="_")
  temp$POP<-pops[p]
  temp$POP_COLOR<-cols[p]
  names(temp)[2]<-"CR"
  temp<-temp[order(-temp$CR),]
  gcr.all.dat<-rbind(gcr.all.dat, temp)
}

gcr.all.dat10<-subset(gcr.all.dat, GENE=="HELLO")
for(p in 1:length(pops)){
  temp<-subset(gcr.all.dat, POP==pops[p])
  temp<-temp[order(-temp$CR),]
  temp<-temp[c(1:10),]
  gcr.all.dat10<-rbind(gcr.all.dat10, temp)
}

pdf("GCR_bypop_top10.dotplot.pdf", width=21, useDingbats=FALSE)
plot(x=seq(1,nrow(gcr.all.dat10)), y=gcr.all.dat10$CR, pch=16, cex=2, cex.lab=1.8, cex.axis=1.4, col=gcr.all.dat10$POP_COLOR, xaxt="n", xlab="", ylab="Gene Carrier Rate (GCR)", ylim=c(0, round(max(dat_01$MAX),2)))
axis(side=1, at=seq(1,nrow(gcr.all.dat10)), labels=gcr.all.dat10$POP_GENE, las=2, cex.axis=1.4)
for(i in seq(0.01, 0.12, 0.01)){
  abline(h=i,  lty=3)
}
dev.off()



#Plot pan-ethnic gene frequencies for top 30 genes
#Figure S1
dat_top30<-dat[order(-dat$MAX_AF),]
dat_top30<-dat_top30[c(1:30),]
pdf("GCR_allpop_top30.dotplot.pdf", width=10, height=10, useDingbats=FALSE)
plot(x=seq(1,nrow(dat_top30)), y=dat_top30$AFR, pch=16, cex=2, cex.axis=1.4, cex.lab=1.8, col=cols[1], xaxt="n", xlab="", ylab="Gene Carrier Rate (GCR)",, ylim=c(0, round(max(dat_top30$MAX),2)))
axis(side=1, at=seq(1,nrow(dat_top30)), labels=dat_top30$GENE, las=2)
points(x=seq(1,nrow(dat_top30)), y=dat_top30$AMR, cex=2, pch=16, col=cols[2])
points(x=seq(1,nrow(dat_top30)), y=dat_top30$ASJ, cex=2, pch=16, col=cols[3])
points(x=seq(1,nrow(dat_top30)), y=dat_top30$EAS, cex=2, pch=16, col=cols[4])
points(x=seq(1,nrow(dat_top30)), y=dat_top30$NFE, cex=2, pch=16, col=cols[5])
points(x=seq(1,nrow(dat_top30)), y=dat_top30$SAS, cex=2, pch=16, col=cols[6])
abline(h=0.01, lty=2, col="black")
legend(x=nrow(dat_top30)*0.9, y=max(dat_top30$MAX), legend=c("AFR", "AMR", "ASJ", "EAS", "NFE", "SAS"),
       col=cols, pch=16, cex=0.8)
dev.off()



#CCR dot plot
#Figure 3A
ccr.cumm.dat<-data.frame(GENE=character(), GCR=numeric(), CCR=numeric(), POP=character(), stringsAsFactors = F)
for(p in 1:length(pops)){
  temp<-dat[,c(1, which(names(dat)==pops[p]))]
  names(temp)[2]<-"GCR"
  temp<-temp[order(-temp$GCR),]
  temp$CCR<-0
  temp[1,]$CCR<-temp[1,]$GCR
  for(i in 2:nrow(temp)){
   temp[i,]$CCR<-(1-prod(1-temp[c(1:i),]$GCR))
  }
  temp$POP<-pops[p]
  ccr.cumm.dat<-rbind(ccr.cumm.dat, temp)
}

pdf("CCR_specific_cummulative.dotplot.pdf", useDingbats=FALSE)
  plot(1, type="n", xlab="Number of Genes Screened", ylab="Cummulative Carrier Rate (CCR)", xlim=c(0, nrow(ccr.cumm.dat[ccr.cumm.dat$POP=="AFR",])), ylim=c(0, round(max(ccr.cumm.dat$CCR)+0.05,1)), cex.axis=1.4, cex.lab=1.8)
  for(i in 1:length(pops)){
    points(x=seq(1,nrow(ccr.cumm.dat[ccr.cumm.dat$POP==pops[i],])), y=ccr.cumm.dat[ccr.cumm.dat$POP==pops[i],]$CCR, pch=16, col=cols[i])
  }
  legend(x=0.9*nrow(ccr.cumm.dat)/length(pops), y=0.2, legend=pops, col=cols, pch=16, cex=0.8)
dev.off()



#Barplot of CCRs
#Figure 3B
ccr.spec.sum.dat<-data.frame(POP=pops, CCR=0,CCR_001=0,CCR_005=0, CCR_01=0,stringsAsFactors = F)
for(i in 1:nrow(ccr.spec.sum.dat)){
  temp<-dat[,c("GENE",ccr.spec.sum.dat[i,]$POP)]
  names(temp)[2]<-c("GCR")
  ccr.spec.sum.dat[i,]$CCR<-(1-prod(1-as.numeric(temp[,2]), na.rm=T))
  ccr.spec.sum.dat[i,]$CCR_001<-(1-prod(1-as.numeric(temp[temp$GCR>0.001,c(2)]), na.rm=T))
  ccr.spec.sum.dat[i,]$CCR_005<-(1-prod(1-as.numeric(temp[temp$GCR>0.005,c(2)]), na.rm=T))
  ccr.spec.sum.dat[i,]$CCR_01<-(1-prod(1-as.numeric(temp[temp$GCR>0.01,c(2)]), na.rm=T))
}

plot.ccr.spec.sum.dat<-melt(ccr.spec.sum.dat)
names(plot.ccr.spec.sum.dat)<-c("POP", "FILTER", "CCR")
pdf("CCR_specific_sum.barplot.pdf")
ggplot(plot.ccr.spec.sum.dat, aes(x=as.factor(POP), y=CCR, fill=FILTER)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
 scale_fill_manual(values=cols2)+
  xlab("Ancestry")+ylab("Cummulative Carrier Rate (CCR)")+
  theme_bw()
dev.off()

#Bar plot of number of genees
#Figure 4A 
genes.count.dat<-data.frame(POP=pops,genes001=0,genes005=0, genes01=0,stringsAsFactors = F)
for(i in 1:nrow(genes.count.dat)){
  genes.count.dat[i,]$genes001<-nrow(dat[dat[,which(names(dat)==genes.count.dat[i,]$POP)]>0.001,])
  genes.count.dat[i,]$genes005<-nrow(dat[dat[,which(names(dat)==genes.count.dat[i,]$POP)]>0.005,])
  genes.count.dat[i,]$genes01<-nrow(dat[dat[,which(names(dat)==genes.count.dat[i,]$POP)]>0.01,])
}
genes.count.dat[nrow(genes.count.dat),]<- c("X_PAN" ,nrow(dat_001), nrow(dat_005), nrow(dat_01))

genes.count.plot.dat<-melt(genes.count.dat, id.vars="POP")
names(genes.count.plot.dat)<-c("POP", "FILTER", "GENES")
genes.count.plot.dat$GENES<-as.numeric(genes.count.plot.dat$GENES)
genes.count.plot.dat$POP<-as.factor(genes.count.plot.dat$POP)
pdf("gene_counts_by_ancestry.barplot.pdf")
ggplot(genes.count.plot.dat, aes(x=POP, y=GENES, fill=FILTER)) + ylab("Number of Genes")+xlab("Ancestry")+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  scale_fill_manual(values=cols2[2:4])+
  theme_bw()
dev.off()


#Cummulative carrier rate using pan-ethnic panels at various thresholds
#Figure 4B-D: 
ccr.sum.pan.dat<-data.frame(POP=pops, SPECIFIC=0,PAN_ETHNIC=0,stringsAsFactors = F)
thresholds<-c(0.01, 0.005, 0.001)
for(t in 1:length(thresholds)){
treshold<-thresholds[t]
for(i in 1:nrow(ccr.sum.pan.dat)){
  temp<-dat[,c("GENE",ccr.sum.pan.dat[i,]$POP, "MAX_AF"),]
  names(temp)[2]<-c("GCR")
  ccr.sum.pan.dat[i,]$SPECIFIC<-(1-prod(1-as.numeric(temp[temp$GCR>threshold,]$GCR), na.rm=T))
  ccr.sum.pan.dat[i,]$PAN_ETHNIC<-(1-prod(1-as.numeric(temp[temp$MAX_AF>threshold,]$GCR), na.rm=T))
}

ccr.sum.pan.plot.dat<-melt(ccr.sum.pan.dat, id.vars="POP")
names(ccr.sum.pan.plot.dat)<-c("POP", "FILTER", "CCR")
pdf(paste("CCR_panethnic_", threshold, ".barplot.pdf", sep=""))
ggplot(ccr.sum.pan.plot.dat, aes(x=as.factor(POP), y=CCR, fill=FILTER))+ylim(c(0,0.75))+
  geom_bar(position=position_dodge(), stat="identity", colour='black')+
  xlab("Ancestry")+ylab("Cummulative Carrier Rate (CCR)")+
  scale_fill_manual(values=c("darkorange", "blue"))+
  theme_bw()
dev.off()
}


#At risk couple rate
#Figure 4A-D
corr.rec.dat<-dat_001 #repeat using dat, dat_01, dat_005, and dat_001

#Make Empty matrix
products<-matrix(nrow=7,ncol=7)
rownames(products)<-pops
colnames(products)<-pops

#Find interpopulation products
for(i in 1:nrow(products)){
  for(j in 1:ncol(products)){
    pop1=rownames(products)[i]
    pop2=rownames(products)[j]
    products[i,j]<-(1-prod(1-(as.vector(corr.rec.dat[,c(pop1)])*as.vector(corr.rec.dat[,c(pop2)]))))
  }
}

products<-round(products*10000, 0)
library(corrplot)
pdf("couple_carrier_rate_001_genes.pdf")
corrplot(products, method="color", type="upper", is.corr=F,addCoef.col = "black",diag=T,cl.lim=c(0,300),tl.col="black",addgrid.col = "black")
dev.off()
