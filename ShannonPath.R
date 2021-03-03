files<-list.files("/Users/hamiji/Downloads/untitled folder", full.names = T)

all_gt<-matrix(0, nrow = 5, ncol = length(files))
for(i in 1:length(files)){
  vcf<-read.vcfR(files[i])
  gt<-extract.gt(vcf)
}

colnames(all_gt)<-all_gt[1,]
all_gt<-all_gt[-1,]
rownames(all_gt)<-c("maltose", "palatinose", "galactose")
write.csv(all_gt, "SNPs.csv")


for(i in 1:length(files)){
  vcf<-read.vcfR(files[i])
  gt<-extract.gt(vcf)
  maltose<-c(length(which(grepl("Cluster5_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster1_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster0_", rownames(gt))=="TRUE")))
  palatinose<-c(length(which(grepl("Cluster7_", rownames(gt))=="TRUE")),
                length(which(grepl("Cluster9_", rownames(gt))=="TRUE")),
                length(which(grepl("Cluster11_", rownames(gt))=="TRUE")),
                length(which(grepl("Cluster17_", rownames(gt))=="TRUE")),
                length(which(grepl("Cluster20_", rownames(gt))=="TRUE")),
                length(which(grepl("Cluster21_", rownames(gt))=="TRUE")),
                length(which(grepl("Cluster28_", rownames(gt))=="TRUE")))
galactose<-c(length(which(grepl("Cluster3_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster4_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster2_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster24_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster25_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster29_", rownames(gt))=="TRUE")),
                    length(which(grepl("Cluster12_", rownames(gt))=="TRUE")))
maltotriose<-c(length(which(grepl("Cluster5_", rownames(gt))=="TRUE")),
               length(which(grepl("Cluster1_", rownames(gt))=="TRUE")),
               length(which(grepl("Cluster26_", rownames(gt))=="TRUE")))
maltose<-maltose/sum(maltose)
maltose[which(maltose == 0)]<-1
maltose_H<--sum(maltose*log2(maltose))
palatinose<-palatinose/sum(palatinose)
palatinose[which(palatinose == 0)]<-1
palatinose_H<--sum(palatinose*log2(palatinose))
galactose<-galactose/sum(galactose)
galactose[which(galactose == 0)]<-1
galactose_H<--sum(galactose*log2(galactose))
maltotriose<-maltotriose/sum(maltotriose)
maltotriose[which(maltotriose == 0)]<-1
maltotriose_H<--sum(maltotriose*log2(maltotriose))
sample<-strsplit(colnames(gt), "/")[[1]][1]
all_gt[1:5,i]<-c(sample, maltose_H, palatinose_H, galactose_H, maltotriose_H)
}

colnames(all_gt)<-all_gt[1,]
all_gt<-all_gt[-1,]
rownames(all_gt)<-c("maltose", "palatinose", "galactose", "maltotriose")
write.csv(all_gt, "SNPs_H.csv")
