## DEseq experiments: 
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
## to write table of TPM counts
table_counts<-counts(ds_all)
write.csv(table_counts,"~/posdoc/florida_trans_spatial_analyses-dir/supplement_final/tables/table_S3_unigene_counts.csv",sep = "\t")
head(table_counts)
#### Mean by sample type: 
table_means<-as.data.frame(table_counts[,c(1:3)])
table_means$`mean_WW-7-M`<-as.data.frame(rowMeans(table_counts[,grep('WW-7-M',colnames(table_counts))]))[,1]
table_means$`mean_WW-7-BS`<-as.data.frame(rowMeans(table_counts[,grep('WW-7-BS',colnames(table_counts))]))[,1]
table_means$`mean_WW-23-M`<-as.data.frame(rowMeans(table_counts[,grep('WW-23-M',colnames(table_counts))]))[,1]
table_means$`mean_WW-23-BS`<-as.data.frame(rowMeans(table_counts[,grep('WW-23-BS',colnames(table_counts))]))[,1]
table_means$`mean_D-7-M`<-as.data.frame(rowMeans(table_counts[,grep('D-7-M',colnames(table_counts))]))[,1]
table_means$`mean_D-7-BS`<-as.data.frame(rowMeans(table_counts[,grep('D-7-BS',colnames(table_counts))]))[,1]
table_means$`mean_D-23-M`<-as.data.frame(rowMeans(table_counts[,grep('D-23-M',colnames(table_counts))]))[,1]
table_means$`mean_D-23-BS`<-as.data.frame(rowMeans(table_counts[,grep('D-23-BS',colnames(table_counts))]))[,1]
table_means<-table_means[,-c(1:3)]
#####################
> resultsNames(ds_all)
[1] "Intercept"          "batch_B_vs_A"       "batch_C_vs_A"       "plant_PO2_vs_PO1"  
[5] "time_7_vs_23"       "tissue_M_vs_BS"     "treatment_D_vs_WW"  "tissueM.treatmentD"

########
[1] "Intercept" 
"tissue_M_vs_BS"    (DE M vs BS) 
"treatment_D_vs_WW" (DE GENES FROM WW TO D) 
"tisseM.treatmentD" (DE GENES FROM WW TO D in the messophyl) 
### WW is the reference level for treatment
### BS is the reference level for tissue
### 23 is the reference level for time
res_tissue <- as.data.frame(results(ds_all, contrast=c("tissue","M","BS")))
colnames(res_tissue)<-paste("MvsBS", colnames(res_tissue), sep = "_")
res_time = as.data.frame(results(ds_all, contrast=c("time","7","23")))
colnames(res_time)<-paste("7vs23", colnames(res_time), sep = "_")
res_treatment = as.data.frame(results(ds_all, contrast=c("treatment","D","WW")))
colnames(res_treatment)<-paste("DvsWW", colnames(res_treatment), sep = "_")
res_plant = as.data.frame(results(ds_all, contrast=c("plant","PO1","PO2")))
colnames(res_plant)<-paste("PO1vsPO2", colnames(res_plant), sep = "_")

deseq_all_factors_tab<-cbind(res_tissue,res_treatment,res_time,res_plant)

###########################
## get stat of each contrast: 

nrow(deseq_all_factors_tab)
## M vs BS

nrow(deseq_all_factors_tab[ which(deseq_all_factors_tab$MvsBS_padj < 0.05), ])
MvBS_sig<-deseq_all_factors_tab[ which(deseq_all_factors_tab$MvsBS_padj < 0.05), ]
#n genes up in BS
nrow(MvBS_sig[which(MvBS_sig$MvsBS_log2FoldChange < 0),])
1030*100/2702
#n genes up in BS
nrow(MvBS_sig[which(MvBS_sig$MvsBS_log2FoldChange > 0),])

## WW vs D
nrow(deseq_all_factors_tab[ which(deseq_all_factors_tab$DvsWW_padj < 0.05), ])
1196*100/22509
DvsWW_sig<-deseq_all_factors_tab[ which(deseq_all_factors_tab$DvsWW_padj < 0.05), ]
#n genes up in WW
nrow(DvsWW_sig[which(DvsWW_sig$DvsWW_log2FoldChange < 0),])
451*100/1196
#n genes up in BS
nrow(DvsWW_sig[which(DvsWW_sig$DvsWW_log2FoldChange > 0),])

## 23 vs 7
nrow(deseq_all_factors_tab[ which(deseq_all_factors_tab$`7vs23_padj` < 0.05), ])
7513*100/22509
h7vs23_sig<-deseq_all_factors_tab[ which(deseq_all_factors_tab$`7vs23_padj` < 0.05), ]
#n genes up in 23
nrow(h7vs23_sig[which(h7vs23_sig$`7vs23_log2FoldChange` < 0),])
3553*100/7513
#n genes up in 7
nrow(h7vs23_sig[which(h7vs23_sig$`7vs23_log2FoldChange` > 0),])

## PLANTS

## 23 vs 7
nrow(deseq_all_factors_tab[ which(deseq_all_factors_tab$PO1vsPO2_padj < 0.05), ])
862*100/22509


#nrow(deseq_all_factors_tab[!duplicated(deseq_all_factors_tab[,'unigenes']),])
nrow(deseq_all_factors_tab[ which(deseq_all_factors_tab$MvsBS_padj < 0.05), ])


sig_M_effect_to_D_day <- sig_M_effect_to_D_day [ order(sig_M_effect_to_D_day$padj), ]
summary(sig_M_effect_to_D_day)


## CCM  
sig_M_effect_to_D_day<-as.data.frame(sig_M_effect_to_D_day)
sig_M_effect_to_D_day$unigenes<-rownames(sig_M_effect_to_D_day)
sig_M_effect_to_D_annot<-merge(sig_M_effect_to_D_day,comb_annot, by ='unigenes', all = FALSE)
sig_M_effect_to_D_annot<-sig_M_effect_to_D_annot[!duplicated(sig_M_effect_to_D_annot[,'unigenes']),]
View(sig_M_effect_to_D_annot)
#### THAT WAS DIFFERENCES BETWEEN TISSUES

### NOW, genes affected 
res_treatment <- results(ds_all, contrast=c("treatment","D","WW"))
#res_treatment = results(ds_all, list( c("treatment_D_vs_WW") ))
sig_treatment <- res_treatment[ which(res_treatment$padj < 0.05), ]
sig_treatment <- sig_treatment [ order(sig_treatment$padj), ]
summary(sig_treatment)
## CCM  
sig_treatment<-as.data.frame(sig_treatment)
sig_treatment$unigenes<-rownames(sig_treatment)
sig_treatment_annot<-merge(sig_treatment,comb_annot, by ='unigenes', all = FALSE)
sig_treatment_annot<-sig_treatment_annot[!duplicated(sig_treatment_annot[,'unigenes']),]
View(sig_treatment_annot)


