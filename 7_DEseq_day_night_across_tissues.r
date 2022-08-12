### DESEQ DAY AND NIGHT: 
comb_annot<-read.csv("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/tableX_phyloannotation_and_Ian_annotations.csv",sep = "\t")
comb_annot[,1]<-NULL
so_norm_filt_matrix<-read.csv("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/sleuth_gene_mode_13may21/kallisto_sleuth_tpm_matrix_for_SI.csv", row.names = 1)
colnames(so_norm_filt_matrix)<-gsub('X','',colnames(so_norm_filt_matrix))
colnames(so_norm_filt_matrix)<-gsub('[.]','-',colnames(so_norm_filt_matrix))
counts3<-round(so_norm_filt_matrix, digits = 0)
#### samples info
library(stringr)
g<-(str_split_fixed(as.list(colnames(counts3)), "_",2)[,2])
g<-gsub("DD", "D", g)
samples <-data.frame(treatment = str_split_fixed(g, "-",3)[,2])
samples$tissue<-as.factor(str_split_fixed(g, "-",4)[,4])
samples$time<-as.factor(str_split_fixed(g, "-",4)[,3])
samples$plant<-as.factor(str_split_fixed(g, "-",4)[,1])
samples$library<-colnames(counts3)
samples$library<-gsub('-','-',samples$library)
samples$tissuetimetreat <-paste(samples$tissue,samples$time,samples$treatment, sep = '_' )
batch<-(str_split_fixed((str_split_fixed(as.list(colnames(counts3)), "_",2)[,1]),"-",2)[,2])
samples$batch<-as.factor(gsub("^.*-", "", batch))
samples$tissuetime <-paste(samples$tissue,samples$time, sep = '_' )

samples$treatment <- as.factor(samples$treatment)
samples$treatment <- relevel(samples$treatment, "WW")
#####################
library(DESeq2)
ds <- DESeqDataSetFromMatrix(countData=counts3, colData=samples, design=~ batch + plant + time + tissue + treatment)
ds$treatment <- relevel(samples$treatment, "WW")
#ds$tissuetimetreat <- relevel(ds$tissuetimetreat, ref = "BS_23_WW")
#ds$tissuetime <- relevel(ds$tissuetime, ref = "BS_23")
keep <- rowSums(counts(ds)) >= 20
ds <- ds[keep,]
#colnames(ds) <- colnames(counts3)
ds_all <- DESeq(ds)
resultsNames(ds_all)
###
res_tissue <- as.data.frame(results(ds_all, contrast=c("tissue","M","BS")))
colnames(res_tissue)<-paste("MvsBS", colnames(res_tissue), sep = "_")
res_time = as.data.frame(results(ds_all, contrast=c("time","7","23")))
colnames(res_time)<-paste("7vs23", colnames(res_time), sep = "_")
res_treatment = as.data.frame(results(ds_all, contrast=c("treatment","D","WW")))
colnames(res_treatment)<-paste("DvsWW", colnames(res_treatment), sep = "_")
res_plant = as.data.frame(results(ds_all, contrast=c("plant","PO1","PO2")))
colnames(res_plant)<-paste("PO1vsPO2", colnames(res_plant), sep = "_")

deseq_all_factors_tab<-cbind(res_tissue,res_treatment,res_time,res_plant)
######## Night
night_counts<-counts3[,grepl("-23-", colnames(counts3))]
g<-(str_split_fixed(as.list(colnames(night_counts)), "_",2)[,2])
g<-gsub("DD", "D", g)
samples <-data.frame(treatment = str_split_fixed(g, "-",3)[,2])
samples$tissue<-as.factor(str_split_fixed(g, "-",4)[,4])
samples$time<-as.factor(str_split_fixed(g, "-",4)[,3])
samples$plant<-as.factor(str_split_fixed(g, "-",4)[,1])
samples$library<-colnames(night_counts)
samples$library<-gsub('-','-',samples$library)
samples$tissuetimetreat <-paste(samples$tissue,samples$time,samples$treatment, sep = '-' )
batch<-(str_split_fixed((str_split_fixed(as.list(colnames(night_counts)), "_",2)[,1]),"-",2)[,2])
samples$batch<-as.factor(gsub("^.*-", "", batch))

samples$treatment <- as.factor(samples$treatment)
samples$treatment <- relevel(samples$treatment, "WW")
library(DESeq2)
ds_night <- DESeqDataSetFromMatrix(countData=night_counts, colData=samples, design=~plant + tissue + treatment + tissue:treatment)
ds_night$treatment <- relevel(ds_night$treatment, ref = "WW")
#keep <- rowSums(counts(ds_night)) >= 20
ds_night <- ds_night[keep,]
colnames(ds_night) <- colnames(night_counts)
ds_night_drought <- DESeq(ds_night)
res_night_D = results(ds_night_drought, list( c("treatment_D_vs_WW") ))
res_night_BS_effect_to_D = results(ds_night_drought, contrast=c("treatment","D","WW"))
res_night_M_effect_to_D <- results(ds_night_drought, list( c("treatment_D_vs_WW","tissueM.treatmentD") ))
### Day
day_counts<-counts3[,grepl("-7-", colnames(counts3))]
g<-(str_split_fixed(as.list(colnames(day_counts)), "_",2)[,2])
g<-gsub("DD", "D", g)
samples <-data.frame(treatment = str_split_fixed(g, "-",3)[,2])
samples$tissue<-as.factor(str_split_fixed(g, "-",4)[,4])
samples$time<-as.factor(str_split_fixed(g, "-",4)[,3])
samples$plant<-as.factor(str_split_fixed(g, "-",4)[,1])
samples$library<-colnames(day_counts)
samples$library<-gsub('-','-',samples$library)
samples$tissuetimetreat <-paste(samples$tissue,samples$time,samples$treatment, sep = '-' )
batch<-(str_split_fixed((str_split_fixed(as.list(colnames(day_counts)), "_",2)[,1]),"-",2)[,2])
samples$batch<-as.factor(gsub("^.*-", "", batch))

samples$treatment <- as.factor(samples$treatment)
samples$treatment <- relevel(samples$treatment, "WW")
## make gene map for DEseq
library(DESeq2)
ds_day <- DESeqDataSetFromMatrix(countData=day_counts, colData=samples, design=~ plant + tissue + treatment + tissue:treatment)
ds_day$treatment <- relevel(ds_day$treatment, ref = "WW")
#keep <- rowSums(counts(ds)) >= 20
## KEEP COMES FROM DAY AND NIGHT ANALYSIS
ds_day <- ds_day[keep,]
colnames(ds_day) <- colnames(day_counts)
ds_day_drought <- DESeq(ds_day)
######The effect of treatment in BS (the main effect).
res_BS_effect_to_D = results(ds_day_drought, contrast=c("treatment","D","WW"))
######################## The effect of D in Mesophyl
res_M_effect_to_D <- results(ds_day_drought, list( c("treatment_D_vs_WW","tissueM.treatmentD") ))
######### Mesophyll day and night effect: 
meso_counts<-counts3[,grepl("-M", colnames(counts3))]
g<-(str_split_fixed(as.list(colnames(meso_counts)), "_",2)[,2])
g<-gsub("DD", "D", g)
samples <-data.frame(treatment = str_split_fixed(g, "-",3)[,2])
samples$tissue<-as.factor(str_split_fixed(g, "-",4)[,4])
samples$time<-as.factor(str_split_fixed(g, "-",4)[,3])
samples$plant<-as.factor(str_split_fixed(g, "-",4)[,1])
samples$library<-colnames(meso_counts)
samples$library<-gsub('-','-',samples$library)
samples$tissuetimetreat <-paste(samples$tissue,samples$time,samples$treatment, sep = '-' )
batch<-(str_split_fixed((str_split_fixed(as.list(colnames(meso_counts)), "_",2)[,1]),"-",2)[,2])
samples$batch<-as.factor(gsub("^.*-", "", batch))

samples$treatment <- as.factor(samples$treatment)
samples$treatment <- relevel(samples$treatment, "WW")
## make gene map for DEseq
library(DESeq2)
ds_meso <- DESeqDataSetFromMatrix(countData=meso_counts, colData=samples, design=~ plant + treatment + time + treatment:time)
ds_meso$treatment <- relevel(ds_meso$treatment, ref = "WW")
#keep <- rowSums(counts(ds)) >= 20
## KEEP COMES FROM meso AND NIGHT ANALYSIS
ds_meso <- ds_meso[keep,]
colnames(ds_meso) <- colnames(meso_counts)
ds_meso_drought <- DESeq(ds_meso)
######The effect of day in well watered samples (the main effect).
res_ww_effect_to_day = results(ds_meso_drought, contrast=c("time","7","23"))
######################## The effect of day in drought samples
res_drought_effect_to_day <- results(ds_meso_drought, list( c("time_7_vs_23","treatmentD.time7") ))
