dat<-read.delim("clinvar_subset_table.txt", header=T, stringsAsFactors = F, sep="\t")
dat$SNP<-paste(dat$chrom, dat$pos, dat$ref, dat$alt, sep=":")
dat$AF_POPMAX<-as.numeric(dat$AF_POPMAX)

#Subset to severe recessive genes
severe<-read.delim("wasserman_gene_list.txt", header=F, stringsAsFactors = F, sep="\t")
dat<-subset(dat, symbol%in%severe$V1 & symbol!="ABCC6")

#Calculate carrier rates for each variant 
pops<-c("ALL", "AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "SAS")
dat$ALL_CR<-(dat$AC-dat$Hom)/(dat$AN/2)
for( i in 2:length(pops)){
  dat$newcol<-(dat[,paste("AC_", pops[i], sep="")]- dat[,paste("Hom_", pops[i], sep="")])/(dat[,paste("AN_", pops[i], sep="")]/2)
  names(dat)[ncol(dat)]<-paste(pops[i], "_CR", sep="")
}

#Calculate USA composite carrier rates
dat$USA_CR<-0.678*dat$NFE_CR+0.019*dat$ASJ_CR+0.122*dat$AFR_CR+0.0175*dat$EAS_CR+0.00997*dat$SAS_CR+0.1533*dat$AMR_CR

#remove bad snps
faulty<-read.delim("Table_S1.txt", header=T, stringsAsFactors=F, sep="\t")
faulty$SNP<-paste(faulty$chrom, faulty$pos, faulty$ref, faulty$alt, sep=":")

dat$SNP<-paste(dat$chrom, dat$pos, dat$ref, dat$alt, sep=":")
dat<-dat[ !(dat$SNP%in%faulty$SNP),]

#Calculate GCRs for each gene
sum.dat<-data.frame(symbol=unique(dat$symbol), ALL=0, AFR=0,AMR=0, ASJ=0, EAS=0, FIN=0, NFE=0, SAS=0, USA=0, stringsAsFactors = F)
for(i in 1:nrow(sum.dat)){
  for(p in 2:ncol(sum.dat)){
    pop<-names(sum.dat)[p]
    sum.dat[i,pop]<-(1-prod(1-as.numeric(dat[dat$symbol==sum.dat[i,]$symbol,paste(pop,"_CR", sep="")]), na.rm=T))
  }
}
write.table(sum.dat, "clinvar_subset_summed_table.txt", row.names=F, col.names=T, sep="\t", quote=F)
