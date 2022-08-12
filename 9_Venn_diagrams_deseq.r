## Venn diagram ###
library("ggVennDiagram")
head(deseq_all_factors_tab)
## get significant genes per variable
library(dplyr)
plantPO2_sig <- dplyr::filter(deseq_all_factors_tab, PO1vsPO2_padj < 0.05)
plantPO2_sig_up <- dplyr::filter(plantPO2_sig, PO1vsPO2_log2FoldChange > 0)
plantPO2_sig_down <- dplyr::filter(plantPO2_sig, PO1vsPO2_log2FoldChange < 0)

treatmentD_sig <- dplyr::filter(deseq_all_factors_tab, DvsWW_padj <= 0.05)
treatmentD_sig_up <- dplyr::filter(treatmentD_sig, DvsWW_log2FoldChange > 0)
treatmentD_sig_down <- dplyr::filter(treatmentD_sig, DvsWW_log2FoldChange < 0)

tisseM_sig <- dplyr::filter(deseq_all_factors_tab, MvsBS_padj <= 0.05)
tisseM_sig_up <- dplyr::filter(tisseM_sig, MvsBS_log2FoldChange > 0)
tisseM_sig_down <- dplyr::filter(tisseM_sig, MvsBS_log2FoldChange < 0)

time7_sig <- dplyr::filter(deseq_all_factors_tab, `7vs23_padj` <= 0.05)
time7_sig_up <- dplyr::filter(time7_sig, `7vs23_log2FoldChange` > 0)
time7_sig_down <- dplyr::filter(time7_sig, `7vs23_log2FoldChange` < 0)


## CALCULATE STATISTICS) 
(nrow(plantPO2_sig)/nrow(deseq_all_factors_tab))*100
(nrow(time7_sig_down)/nrow(time7_sig))*100

##### VENN DIAGRAM


x <- list(
  Treatment = treatmentD_sig$unigenes,
  Tissue = tisseM_sig$unigenes,
  Time = time7_sig$unigenes,
  Individual = plantPO2_sig$unigenes
)



up <- list(
  Drought = treatmentD_sig_up$unigenes,
  Mesohpyll = tisseM_sig_up$unigenes,
  Day = time7_sig_up$unigenes,
  Plant2 = plantPO2_sig_up$unigenes
)


down <- list(
  Watered = treatmentD_sig_down$unigenes,
  Bundle_Sheath = tisseM_sig_down$unigenes,
  Night = time7_sig_down$unigenes,
  Plant = plantPO2_sig_down$unigenes
)
pdf('~/Desktop/figures_lcm_final/venn_diagrams_version2.pdf',width = 6, height = 6)

ggVennDiagram(x, label_alpha = 0)
ggVennDiagram(up, label_alpha = 0)
ggVennDiagram(down, label_alpha = 0)
dev.off()
#### extract intersect of CCM genes: 
Reduce(intersect,  list(treatmentD_sig_up$unigenes,
tisseM_sig_up$unigenes,time7_sig_up$unigenes))
DE_all_groups<-Reduce(intersect, x[-4])


########## STATISTICS NUMBER OF GENES DE


## % reads mapped to CCM genes#### 
kall_read_abundances<-kallisto_table(so, use_filtered = FALSE, normalized = FALSE,
                                     include_covariates = TRUE)
head(kall_read_abundances)
library(tidyr)
### long to wide (matrix of counts)
raw_read_matrix <- spread(kall_read_abundances[,c(1,2,6)], sample, est_counts)
## add unigenes naming
raw_read_matrix$unigenes = sub("_i.*", "", as.character(raw_read_matrix$target_id))
## (ccm reads in M libraries/total reads all genes in M libraries)*100
total_7M<-colSums(raw_read_matrix[,grep('7-M',names(raw_read_matrix))])
#### subset to CCM
ccm_raw_reads<-raw_read_matrix[raw_read_matrix$unigenes %in% comb_annot$unigenes,]
total_ccm_7M<-colSums(ccm_raw_reads[,grep('7-M',names(ccm_raw_reads))])
################
mean((total_ccm_7M/total_7M)*100)
sd((total_ccm_7M/total_7M)*100)
## 10.97285 sd 2.488944
##################### BS

total_7BS<-colSums(raw_read_matrix[,grep('D-7|W-7',names(raw_read_matrix))])
#### subset to CCM
total_ccm_7BS<-colSums(ccm_raw_reads[,grep('D-7|W-7',names(ccm_raw_reads))])
mean((total_ccm_7BS/total_7BS)*100)
sd((total_ccm_7BS/total_7BS)*100)
12.18361 % 2.938403

