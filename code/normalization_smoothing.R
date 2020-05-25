# clear data
rm(list=ls())

## Z score normalization
# Go to the directory containing data files and import them into R by running the following commands:
setwd("/Users/luliu/Box/Bioinfo_analysis/FUS_Replication timing/Replication timing_08312018")
getwd()
chrs<-append(paste0("chr",rep(1:22)),"chrX")
for (sample in c("Clone110_RT_R1-X","Clone110_RT_R2-X","FUSClone110_RT_R1-X","FUSClone110_RT_R2-X","U2OS_RT_R1-X","U2OS_RT_R2-X")){
  ## Z score nomorlization
    # read data
    # sample<-"Clone110_RT_R1-X"
    sample_data<-read.table(paste0("./data/",sample,".bedgraph",sep="") ,header = FALSE)
    str(sample_data)
    # remove other chr, only keep the chr1 to chr22 and chrX
    sample_data<-sample_data[sample_data$V1 %in% chrs,]
    # z score nomorlization
    sample_data[,4]<- (sample_data[,4]-mean(sample_data[,4]))/sd(sample_data[,4])
    summary(sample_data)
    # explore table
    write.table(sample_data, paste0("./01_02_nomorlization_smoothing/Z_Score_normalization/",sample,"_z_score_normal.bedgraph",sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
    
  ## Loess smoothing
    # Initialize an R-list to stock data sets
    AllLoess=list()
    #Perform Loess smoothing
    for(i in 1:(ncol(sample_data)-3)){
      AllLoess[[i]]=data.frame();
      cat("Current dataset:", colnames(sample_data)[i+3], "\n");
      for(Chr in chrs){
        RTb=subset(sample_data, sample_data$V1==Chr);
        lspan=2000000/(max(RTb$V2)-min(RTb$V2));
        cat("Current chrom:" , Chr, "\n");
        RTla=loess(RTb[,i+3] ~ RTb$V2, span=lspan);
        RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x, sample_data[which( sample_data$V1==Chr & sample_data$V2 %in% RTla$x),3],RTla$fitted);
        colnames(RTl)=c("chr" , "start" , "end" ,"counts_norm");
        if(length(AllLoess[[i]])!=0){
          AllLoess[[i]]=rbind(AllLoess[[i]],RTl)};
        if(length(AllLoess[[i]])==0){
          AllLoess[[i]] = RTl}
       }
    }

    ## Register the Loess-smoothed data into bedGraph files
  for(i in 1:length(AllLoess)){
    write.table(AllLoess[[i]][complete.cases(AllLoess[[i]]),], gsub(".bg" , "Loess.bedGraph" , paste0("./01_02_nomorlization_smoothing/smoothing/",sample,"_Loess_smoothing.bedGraph",sep="")), sep= "\t" , row.names=FALSE, quote=FALSE, col.names = FALSE)
    }
}   
    










