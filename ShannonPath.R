files<-list.files("/Users/hamiji/Downloads/Shannon_Index", full.names = T)
files_names<-list.files("/Users/hamiji/Downloads/Shannon_Index")


all_gt<-data.frame()


GS_clusters<-read.csv("/Users/hamiji/GS_Clusters.csv", row.names = 1)


for(i in 1:length(files)){
  vcf<-read.vcfR(files[i])
  gt<-extract.gt(vcf, return.alleles = F)
  gt_df<-do.call(rbind.data.frame, strsplit(gt, "/"))
  rownames(gt_df)<-rownames(gt)
  
  SNPs_cnt<-data.frame()
  for(j in 1:nrow(gt_df)){
    SNPs<-length(unique(unlist(gt_df[j,])))
    SNPs_cnt<-rbind(SNPs_cnt,SNPs)
  }
  
  gt_df<-cbind(gt_df, SNPs_cnt)
  
  SNP_go_total<-data.frame()
  for(j in 1:nrow(GS_clusters)){
    GO_subset<-paste0(unique(unlist(GS_clusters[j,])),"_")[-length(paste0(unique(unlist(GS_clusters[j,])),"_"))]
    SNP_sum<-cbind(unlist(lapply(GO_subset, function(x){
      sum(gt_df[which(grepl(x, rownames(gt_df))=="TRUE"),ncol(gt_df)])
    })))
    GO_H<-SNP_sum/sum(SNP_sum)
    GO_H[which(GO_H==0)]<-1
    GO_H_H<--sum(GO_H*log2(GO_H))
    SNP_go_total<-rbind(SNP_go_total,GO_H_H)
  }
  
  
  y<-apply(GS_clusters, 1, function(x){
    GO_subset<-paste0(unique(unlist(x)),"_")[-length(paste0(unique(unlist(x)),"_"))]
    SNP_sum<-cbind(unlist(lapply(GO_subset, function(x){
      sum(gt_df[which(grepl(x, rownames(gt_df))=="TRUE"),ncol(gt_df)])
    })))
    GO_H<-SNP_sum/sum(SNP_sum)
    GO_H[which(GO_H==0)]<-1
    GO_H_H<--sum(GO_H*log2(GO_H))
  })
  
  
  
  maltose<-maltose/sum(maltose)
  maltose[which(maltose == 0)]<-1
  maltose_H<--sum(maltose*log2(maltose))
  
  
}

colnames(all_gt)<-files_names
rownames(all_gt)<-c("maltose", "palatinose", "galactose", "maltotriose")
write.csv(all_gt, "/Users/hamiji/SNPs_H.csv")
