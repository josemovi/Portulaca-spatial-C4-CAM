### RECONCILE TABLES: 
# table containing all results desesq using all variables
deseq_all_factors_tab
# tables day night separation: 
res_BS_effect_to_D<-as.data.frame(res_BS_effect_to_D)
# add variable tag 
colnames(res_BS_effect_to_D)<-paste("Day_BS", colnames(res_BS_effect_to_D), sep = "_")
res_M_effect_to_D<-as.data.frame(res_M_effect_to_D)
colnames(res_M_effect_to_D)<-paste("Day_M", colnames(res_M_effect_to_D), sep = "_")
# tables day night separation:
res_night_BS_effect_to_D<-as.data.frame(res_night_BS_effect_to_D)
colnames(res_night_BS_effect_to_D)<-paste("Night_BS", colnames(res_night_BS_effect_to_D), sep = "_")
res_night_M_effect_to_D<-as.data.frame(res_night_M_effect_to_D)
colnames(res_night_M_effect_to_D)<-paste("Night_M", colnames(res_night_M_effect_to_D), sep = "_")

deseq_all_factors_tab<-cbind(deseq_all_factors_tab,res_M_effect_to_D,res_BS_effect_to_D,res_night_M_effect_to_D,res_night_BS_effect_to_D)
#### add mesophyll day and night: 
res_ww_effect_to_day<-as.data.frame(res_ww_effect_to_day)
colnames(res_ww_effect_to_day)<-paste("M_WW_day", colnames(res_ww_effect_to_day), sep = "_")
res_drought_effect_to_day<-as.data.frame(res_drought_effect_to_day)
colnames(res_drought_effect_to_day)<-paste("M_D_day", colnames(res_drought_effect_to_day), sep = "_")
deseq_all_factors_tab<-cbind(deseq_all_factors_tab,res_ww_effect_to_day,res_drought_effect_to_day)


#####################
#######################

deseq_all_factors_tab$unigenes<-rownames(deseq_all_factors_tab)
deseq_all_factors_tab_annot<-merge(deseq_all_factors_tab,comb_annot, by ='unigenes', all = FALSE)
deseq_all_factors_tab_annot<-deseq_all_factors_tab_annot[!duplicated(deseq_all_factors_tab_annot[,'unigenes']),]

##################
##################ANNOTATE TABLE WITH TRINONATE INFO. 
trinotate<-read.delim('~/posdoc/florida_trans_spatial_analyses-dir/trinotate_annotation_report.xls')
trinotate$unigenes = sub("_i.*", "", as.character(trinotate$transcript_id))
trinotate<-trinotate[c('unigenes','sprot_Top_BLASTP_hit')]
### trinotate 
##### 1: get funtional annotations: 
deseq_all_factors_tab_trinot<-merge(deseq_all_factors_tab,trinotate, by= 'unigenes', all=F)
deseq_all_factors_tab_trinot<-deseq_all_factors_tab_trinot[!duplicated(deseq_all_factors_tab_trinot[,1]),]
### reduce columns 
deseq_all_factors_tab_trinot<-deseq_all_factors_tab_trinot[,c(!grepl('lfcSE|stat|pvalue',colnames(deseq_all_factors_tab_trinot)))]
### write table: 
write.csv(deseq_all_factors_tab_trinot,"~/posdoc/florida_trans_spatial_analyses-dir/supplement_final/tables/table_S4_trinotate_deseq.csv",sep = "\t")
### reduce to 500 CCM related genes
deseq_all_factors_tab_trinot_2<-read.csv('~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/table_S4_trinotate_deseq_v1',sep = "\t")
######### merge means
table_means$unigenes<-rownames(table_means)
deseq_all_factors_tab_trinot_2<-merge(deseq_all_factors_tab_trinot_2,table_means, by= 'unigenes')
############
write.csv(deseq_all_factors_tab_trinot_2,"~/posdoc/florida_trans_spatial_analyses-dir/supplement_final/tables/means.csv",sep = "\t")
############  reduce table to the CCM annotated genes
deseq_all_factors_tab_trinot_ccm<-merge(deseq_all_factors_tab_trinot_2,comb_annot, by ='unigenes', all = FALSE)
deseq_all_factors_tab_trinot_ccm<-deseq_all_factors_tab_trinot_ccm[!duplicated(deseq_all_factors_tab_trinot_ccm[,'unigenes']),]
write.csv(deseq_all_factors_tab_trinot_ccm,"~/posdoc/florida_trans_spatial_analyses-dir/supplement_final/tables/Table_S5_CCM_genes.csv",sep = "\t")

deseq_all_factors_tab_trinot_ccm2<-read.csv("~/posdoc/florida_trans_spatial_analyses-dir/supplement_final/tables/Table_S5_CCM_genes_v2.csv",sep = ",")

nrow(deseq_all_factors_tab)
#### Get stats: 
## M vs BS
nrow(deseq_all_factors_tab_trinot_ccm2[ which(deseq_all_factors_tab_trinot_ccm2$MvsBS_padje < 0.05), ])
87*100/170
ccm_MvsBS_sig<-deseq_all_factors_tab_trinot_ccm2[ which(deseq_all_factors_tab_trinot_ccm2$MvsBS_padje < 0.05), ]
#n genes up in BS
nrow(ccm_MvsBS_sig[which(ccm_MvsBS_sig$MvsBS_log2FCe < 0),])
10*100/75
#n genes up in M
nrow(ccm_MvsBS_sig[which(ccm_MvsBS_sig$MvsBS_log2FCe > 0),])


## WW vs D
nrow(deseq_all_factors_tab_trinot_ccm2[ which(deseq_all_factors_tab_trinot_ccm2$DvsWW_padje < 0.05), ])
1196*100/480
ccm_DvsWW_sig<-deseq_all_factors_tab_trinot_ccm2[ which(deseq_all_factors_tab_trinot_ccm2$DvsWW_padje < 0.05), ]
#n genes up in WW
nrow(ccm_DvsWW_sig[which(ccm_DvsWW_sig$DvsWW_log2FCe < 0),])
10*100/75
#n genes up in D
nrow(ccm_DvsWW_sig[which(ccm_DvsWW_sig$DvsWW_log2FCe > 0),])

## 23 vs 7
nrow(deseq_all_factors_tab_trinot_ccm2[ which(deseq_all_factors_tab_trinot_ccm2$X7vs23_padje < 0.05), ])
307*100/480
ccm_h7vs23_sig<-deseq_all_factors_tab_trinot_ccm2[ which(deseq_all_factors_tab_trinot_ccm2$X7vs23_padje < 0.05), ]
#n genes up in 23
nrow(ccm_h7vs23_sig[which(ccm_h7vs23_sig$X7vs23_log2FCe < 0),])
up_night<-ccm_h7vs23_sig[which(ccm_h7vs23_sig$X7vs23_log2FCe < 0),]
View(up_night)
77*100/307
#n genes up in 7
nrow(ccm_h7vs23_sig[which(ccm_h7vs23_sig$X7vs23_log2FCe > 0),])

## PLANTS

## 23 vs 7
nrow(deseq_all_factors_tab[ which(deseq_all_factors_tab$PO1vsPO2_padj < 0.05), ])
862*100/22509


############  reduce table to the CCM annotated genes selected contigs
deseq_all_factors_tab_trinot_ccm_selected<-deseq_all_factors_tab_trinot_ccm2[deseq_all_factors_tab_trinot_ccm2$unigenes %in% deseq_selected_contigs$unigenes,]
write.csv(deseq_all_factors_tab_trinot_ccm_selected,"~/posdoc/florida_trans_spatial_analyses-dir/supplement_final/tables/Table_S6_CCM_selected_contigs.csv",sep = "\t")


#aqui

TRINITY_DN0_c0_g1

deseq_all_factors_tab$unigenes<-rownames(deseq_all_factors_tab)
deseq_all_factors_tab_annot<-merge(deseq_all_factors_tab_trinot,comb_annot, by ='unigenes', all = FALSE)
deseq_all_factors_tab_annot<-deseq_all_factors_tab_annot[!duplicated(deseq_all_factors_tab_annot[,'unigenes']),]




### back to DEseq analysis: 
############# add interesction tissue:treatment
res_dif_response_in_tissue<-as.data.frame(res_dif_response_in_tissue)
colnames(res_dif_response_in_tissue)<-paste("day_tissue_treat", colnames(res_dif_response_in_tissue), sep = "_")
res_dif_response_in_tissue$unigenes<-rownames(res_dif_response_in_tissue)
res_night_dif_response_in_tissue<-as.data.frame(res_night_dif_response_in_tissue)
colnames(res_night_dif_response_in_tissue)<-paste("night_tissue_treat", colnames(res_night_dif_response_in_tissue), sep = "_")
res_night_dif_response_in_tissue$unigenes<-rownames(res_night_dif_response_in_tissue)

deseq_all_factors_tab_annot<-merge(deseq_all_factors_tab_annot,res_dif_response_in_tissue, by ='unigenes', all = FALSE)
deseq_all_factors_tab_annot<-merge(deseq_all_factors_tab_annot,res_night_dif_response_in_tissue, by ='unigenes', all = FALSE)

deseq_all_factors_tab_annot<-deseq_all_factors_tab_annot[!duplicated(deseq_all_factors_tab_annot[,'unigenes']),]

####### RUN ONLY ONCE:: IDENTIFY MODES OF EXPRESSION FOR EACH UNIGENE IN EACH GENE LINEAGE OR 
#### GENE FAMILY AND KEEP ONE REPRESENTATIVE PER EXPRESSION MODE:
# 1) select genes that are DE in any case: 
temp<-deseq_all_factors_tab_annot[,grepl('padj',names(deseq_all_factors_tab_annot))]
temp$row<-rownames(temp)
## select rows when any of the variables are smoller than 0.05
####
library(dplyr)
temp2<-temp %>%  filter_all(any_vars(.< 0.05))
annot_contrast_sig<-deseq_all_factors_tab_annot[temp2$row,]
#########
# 2) order by Pathway, gene family and target_id
annot_contrast_sig<-annot_contrast_sig[with(annot_contrast_sig, order(Pathway, Gene.family,unigenes)),]
nodup_annot_contrast_sig<-annot_contrast_sig[!duplicated(annot_contrast_sig$unigenes),]
nrow(nodup_annot_contrast_sig)
# 3) get only q_values
qval_table<-nodup_annot_contrast_sig[,grepl('padj',names(nodup_annot_contrast_sig))]
rownames(qval_table)<-paste(nodup_annot_contrast_sig$unigenes,nodup_annot_contrast_sig$Gene.family,nodup_annot_contrast_sig$Pamilis.uniprot,nodup_annot_contrast_sig$Athaliana.uniprot,nodup_annot_contrast_sig$phyloannot, sep = '_')
head(qval_table)
# 4) get the beta values
b_table<-nodup_annot_contrast_sig[,grepl('log2',names(nodup_annot_contrast_sig))]
rownames(b_table)<-paste(nodup_annot_contrast_sig$unigenes,nodup_annot_contrast_sig$Gene.family,nodup_annot_contrast_sig$Pamilis.uniprot,nodup_annot_contrast_sig$Athaliana.uniprot,nodup_annot_contrast_sig$phyloannot, sep = '_')
# 5) remove betas that are not significant
####### REMOVE BETAS THAT ARE NOT SIGNIFICANT

##### remove plant o2
b_table3<-b_table[,-4]
###### 
##fit <-kmeans(b_table3, 8)
##### FORCING BETAS INTO GROUPS
b_table4<-b_table3
b_table4[b_table4 < -0.1] <- -1
b_table4[b_table4 > -0.1 & b_table4 < 0.1] <- 0
b_table4[b_table4 > 0.1] <- 1

## remove the betas that are not significant.
b_table4<-replace(b_table4, qval_table[,-4] >0.049,NA) 

########## HOW MANY GROUPS OF EXPRESSION THERE ARE (ignore plant expression) 
b_table4<-as.data.frame(b_table4)
## good order
b_table4<-b_table4[c(2,1,3)]
b_table4$expr_classes<-paste(b_table4$MvsBS_log2FoldChange,b_table4$DvsWW_log2FoldChange,b_table4$`7vs23_log2FoldChange`, sep = '_')

## 6) add expression data
### THE MOS HIHGLY EXPRESSED?
## get contigs names
#library(stringr)
temp<-as.data.frame(str_split_fixed(rownames(b_table4), "_",8)[,c(1,2,3,4,5)])
temp2<-paste(temp$V1,temp$V2,temp$V3,temp$V4,sep = '_')
b_table4$unigenes<-temp2
## subet the tpm matrix to the unigens of interes
###
tpm1<-so_norm_filt_matrix[row.names(so_norm_filt_matrix) %in% temp2,]
##### ADD TPM SUM to each contig: 
csumn<-as.data.frame(colSums(t(tpm1)))
names(csumn)[1]<-'totaltpm'
csumn$unigenes<-row.names(csumn)
### merge losses the rownames, so keep it in column for later
b_table4$unigene_annot<-row.names(b_table4)

b_table5<-merge(b_table4,csumn,by = 'unigenes', sort = FALSE)
b_table4[grepl('NADMDH-6E1',b_table4$unigene_annot),]
qval_table[grepl('NADMDH-6E1',row.names(qval_table)),]

########## HOW MANY GROUPS OF EXPRESSION THERE ARE:
## Filter by expession pattern and by gene family: keep representatives
row.names(b_table5)<- b_table5$unigene_annot
temp<-as.data.frame(str_split_fixed(rownames(b_table5), "_",9)[,c(5:8)])
b_table5$family<-temp$V1
### change those with not family annot but with phylo annotation
#temp<-as.data.frame(str_split_fixed(rownames(b_table5[c(346:359),]), "_",9)[,c(5:8)])
#temp<-as.data.frame(str_split_fixed(temp[,4], "[-]",2))[,1]
#temp
#temp<-gsub("NADMDH","MDH",temp)
#b_table5[c(346:359),'family']<-as.vector(temp)
############ CHECK IF THERE ARE DIFFERENCES IN THE BETAS BETWEEN CONTIGS WITHING EACH FAMILIY
b_table6 <- b_table5[order(b_table5$family, -abs(b_table5$totaltpm) ), ] #sort by id and reverse of abs(value)
## keep the most highly expressed by family and expression mode
b_table6<-b_table6[!duplicated(b_table6[,c('expr_classes','family')]),]
# add family phylo tag
temp<-as.data.frame(str_split_fixed(rownames(b_table6), "_",9)[,c(5:8)])
b_table6$family_phylo<-paste(temp$V1,temp$V4,sep = '_')
###### GET THE HIGHEST BY IN PHYLO, THIS WAY MANY LINEAGES HAVE BEEN DROPPED ALREADY.
b_table7 <- b_table6[order(b_table6$family_phylo, -abs(b_table6$totaltpm) ), ] 
b_table7<-b_table7[!duplicated(b_table7[,'family_phylo']),]
b_table7_annot<-merge(b_table7,comb_annot[,c('unigenes','Pathway')],by ='unigenes', all.x = TRUE, all.y = FALSE)
b_table7_annot<-b_table7_annot[!duplicated(b_table7_annot[,'unigenes']),]
########## add annotations for figure 2 and reorder
##
fig2_anot<-read.csv('~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/pathways_annotation_fig_2.csv', sep = '\t', header = F)
names(fig2_anot)<-c('categories','Pathway')
b_table7_annot2<-merge(b_table7_annot,fig2_anot,by ='Pathway', all.x = TRUE, all.y = FALSE)
b_table7_annot2<-b_table7_annot2[!duplicated(b_table7_annot2[,'unigenes']),]
b_table7_annot2<-b_table7_annot2[with(b_table7_annot2, order(MvsBS_log2FoldChange, DvsWW_log2FoldChange,`7vs23_log2FoldChange`)),]
b_table7_annot2$figure_indexes<-paste(b_table7_annot2$categories,b_table7_annot2$family_phylo, sep = ' ')
write.table(b_table7_annot2,'~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/fig2_DEseq_table_ccm_DE.csv')
##### !!! REMOVED NAs FROM THE FIGURE INDEXES IN OFFICE!!

## CONTINUE HERE AQUI!!!!!!
## all selected contigs:
selected_contigs<-read.csv('~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/fig2_DEseq_table_ccm_manual_reduced_v4.csv')
## reduced selected cointgs for fig2:
selected_contigs<-read.csv('~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/fig2_DEseq_table_ccm_manual_reduced_v8.csv')
selected_contigs<-read.csv('~/Desktop/supplement_final/supplement_barplots_DEseq_table_ccm_manual_reduced.csv', sep = '\t')
head(selected_contigs)
####
# PLOT HEATMAP:
### HERE TO GET DESEQ TABLE REDUCED TO SELECTED CONTIGS
#deseq_selected_contigs<-merge(selected_contigs[,-c(1,3:7)],deseq_all_factors_tab,by = 'unigenes')
deseq_selected_contigs<-merge(selected_contigs[,-c(4:6)],deseq_all_factors_tab,by = 'unigenes')
deseq_selected_contigs<-deseq_selected_contigs[order(match(deseq_selected_contigs$unigenes, selected_contigs$unigenes)), , drop = FALSE]
rownames(deseq_selected_contigs)<-deseq_selected_contigs$indexes
#selected_log<-deseq_selected_contigs[grep('log2',colnames(deseq_selected_contigs))]
#rownames(selected_log)<-deseq_selected_contigs$indexes

#selected_padj<-deseq_selected_contigs[grep('padj',colnames(deseq_selected_contigs))]
#rownames(selected_padj)<-deseq_selected_contigs$indexes

#####


#################### PLOT SIG GENES DIFF RESPONSE M AND B (DAY NIGHT) TO DROUGHT
#library(ggplot2)
#selected_log$genes<-rownames(selected_log)
#selected_log_rev=selected_log[order(nrow(selected_log):1),]
#selected_log_rev$genes <- factor(selected_log_rev$genes, levels = selected_log_rev$genes)
deseq_selected_contigs$genes<-rownames(deseq_selected_contigs)
deseq_selected_contigs_rev=deseq_selected_contigs[order(nrow(deseq_selected_contigs):1),]
deseq_selected_contigs_rev$genes <- factor(deseq_selected_contigs_rev$genes, levels = deseq_selected_contigs_rev$genes)


#pdf("CCM_dif_response_in_tissue_barplot.pdf", width = 10, height = 3)
#
#pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig_supplement_logfc_MvsB.pdf", width = 6, height = 10)
#library(ggplot2)
tissue_p<-ggplot(deseq_selected_contigs_rev, aes(x = genes, y = MvsBS_log2FoldChange, shape = MvsBS_padj < 0.05)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = MvsBS_log2FoldChange > 0 ), show.legend = FALSE) + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(axis.title.y =element_blank()) +
  geom_point(aes(y = 4),
             position = position_dodge(0.9), 
             show.legend = FALSE) +
  scale_shape_manual(values = c(NA, 8)) +
  geom_errorbar(aes(ymin=MvsBS_log2FoldChange-MvsBS_lfcSE, ymax=MvsBS_log2FoldChange+MvsBS_lfcSE), width=.2,
                position=position_dodge(.9)) +
  coord_flip()
#dev.off()


#pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig_supplement_logfc_drought.pdf", width = 6, height = 10)
water_p <- ggplot(deseq_selected_contigs_rev, aes(x = genes, y =  DvsWW_log2FoldChange, shape = DvsWW_padj < 0.05)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = DvsWW_log2FoldChange > 0 ), show.legend = FALSE) + 
#  theme(axis.text.y = element_blank(),axis.title.y =element_blank()) + 
  theme(axis.title.y =element_blank()) +
  geom_point(aes(y = 10),
             position = position_dodge(0.9), 
             show.legend = FALSE) +
  scale_shape_manual(values = c(NA, 8)) +
  geom_errorbar(aes(ymin=DvsWW_log2FoldChange-DvsWW_lfcSE, ymax=DvsWW_log2FoldChange+DvsWW_lfcSE), width=.2,
                position=position_dodge(.9)) +
  coord_flip()
#dev.off()
#################

#pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig_supplement_logfc_time.pdf", width = 6, height = 10)
time_p <-ggplot(deseq_selected_contigs_rev, aes(x = genes, y =  `7vs23_log2FoldChange`, shape = `7vs23_padj` < 0.05)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = `7vs23_log2FoldChange` > 0 ), show.legend = FALSE) + 
#  theme(axis.text.y = element_blank(), axis.title.y =element_blank()) + 
  theme(axis.title.y =element_blank()) +
  geom_point(aes(y = 10),
             position = position_dodge(0.9), 
             show.legend = FALSE) +
  scale_shape_manual(values = c(NA, 8)) +
  geom_errorbar(aes(ymin=`7vs23_log2FoldChange`-`7vs23_lfcSE`, ymax=`7vs23_log2FoldChange`+`7vs23_lfcSE`), width=.2,
                position=position_dodge(.9)) +
  coord_flip()
#dev.off()
#library(cowplot)
p <- plot_grid(tissue_p,time_p,water_p, ncol = 3)
ggsave("~/Desktop/supplement_final/log2_barplots_supplement.pdf", p, width = 12, height = 10)

ggsave("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/Fig_2B_18oct.pdf", p, width = , height = 4)
############################

#dev.off()
ggplot(selected_log_rev, aes(x = genes, y = DvsWW_log2FoldChange)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = DvsWW_log2FoldChange > 0 )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()

ggplot(selected_log_rev, aes(x = genes, y = `7vs23_log2FoldChange`)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = `7vs23_log2FoldChange` > 0 )) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()


plot(selected_log_rev$Day_M_log2FoldChange, 
     selected_log_rev$Night_M_log2FoldChange)
abline(lm(selected_log_rev$Day_M_log2FoldChange ~ selected_log_rev$Night_M_log2FoldChange), col = "red", lwd = 3)
text(paste("Correlation:", round(cor(selected_log_rev$Day_M_log2FoldChange, selected_log_rev$Night_M_log2FoldChange), 2)), x = 25, y = 95)
###### WW samples response to Day Ian style, Day vs Night
#############################
## RESPONSE M TO DROGUTH

#pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig_supplement_logfc_M_to_D.pdf", width = 6, height = 10)
pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig2_C_Mto_time.pdf", width = 4, height = 4.5)

ggplot(deseq_selected_contigs_rev,
       aes(y = genes)) +
  labs(x = "log2FoldChange_M_in_time",y= 'gene') +
  geom_segment(aes(x = deseq_selected_contigs_rev$M_WW_day_log2FoldChange,
                   y = genes,
                   xend = deseq_selected_contigs_rev$M_D_day_log2FoldChange,
                   yend = genes),linetype=2,
               size = 0.5) +
  #  geom_point(aes(x = deseq_selected_contigs_rev$Day_M_log2FoldChange, color = Day_M_log2FoldChange > 0), size = 2, shape = 16) + # DOT
  #  geom_point(aes(x = deseq_selected_contigs_rev$Night_M_log2FoldChange, color = Night_M_log2FoldChange > 0), size = 2, shape = 8) + # STAR
  geom_point(aes(x = deseq_selected_contigs_rev$M_WW_day_log2FoldChange), size = 2, shape = 17, color = 'blue') + # DOT
  geom_point(aes(x = deseq_selected_contigs_rev$M_D_day_log2FoldChange), size = 2, shape = 17, color = 'brown') + # STAR
  geom_point(aes(x = 10, color = M_WW_day_padj <= 0.04999), shape = 8) + # STAR
  geom_point(aes(x = 12, color = M_D_day_padj <= 0.04999), shape = 8) + # STAR
  scale_color_discrete(name = "Transcription") +
  theme_classic() 
#dev.off() 
#aqui

#############################
## RESPONSE M TO DROGUTH

#pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig_supplement_logfc_M_to_D.pdf", width = 6, height = 10)
pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig2_C_Mto_D.pdf", width = 4, height = 4.5)

m_to_d<-ggplot(deseq_selected_contigs_rev,
       aes(y = genes)) +
  labs(x = "log2FoldChange_M_in_D",y= 'gene') +
  geom_segment(aes(x = deseq_selected_contigs_rev$Day_M_log2FoldChange,
                   y = genes,
                   xend = deseq_selected_contigs_rev$Night_M_log2FoldChange,
                   yend = genes),linetype=2,
               size = 0.5) +
#  geom_point(aes(x = deseq_selected_contigs_rev$Day_M_log2FoldChange, color = Day_M_log2FoldChange > 0), size = 2, shape = 16) + # DOT
#  geom_point(aes(x = deseq_selected_contigs_rev$Night_M_log2FoldChange, color = Night_M_log2FoldChange > 0), size = 2, shape = 8) + # STAR
  geom_point(aes(x = deseq_selected_contigs_rev$Day_M_log2FoldChange), size = 3, shape = 17, color = 'blue') + # DOT
  geom_point(aes(x = deseq_selected_contigs_rev$Night_M_log2FoldChange), size = 3, shape = 17, color = 'black') + # STAR
  geom_point(aes(x = 9, color = Day_M_padj <= 0.04999), shape = 8) + # STAR
  geom_point(aes(x = 11, color = Night_M_padj <= 0.04999), shape = 8) + # STAR
      scale_color_discrete(name = "Transcription") +
  theme_classic() 
dev.off()
#theme(legend.position = "bottom")

#####
## RESPONSE BS TO DROGUTH
pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig_supplement_logfc_BS_to_D.pdf", width = 6, height = 10)

b_to_d<-ggplot(deseq_selected_contigs_rev,
       aes(y = genes)) +
  labs(x = "log2FoldChange_B_in_D",y= 'gene') +
  geom_segment(aes(x = deseq_selected_contigs_rev$Day_BS_log2FoldChange,
                   y = genes,
                   xend = deseq_selected_contigs_rev$Night_BS_log2FoldChange,
                   yend = genes),linetype=2,
               size = 0.5) +
  geom_point(aes(x = deseq_selected_contigs_rev$Day_BS_log2FoldChange),size = 3, shape = 17, color = 'blue') + # DOT
  geom_point(aes(x = deseq_selected_contigs_rev$Night_BS_log2FoldChange), size = 3, shape = 17, color = 'black') + # STAR  geom_point(aes(x = 9, color = Day_BS_padj <= 0.04999), shape = 8) + # STAR
  geom_point(aes(x = 9, color = Day_BS_padj <= 0.04999), shape = 8) + # STAR
  geom_point(aes(x = 11, color = Night_BS_padj <= 0.04999), shape = 8) + # STAR
  scale_color_discrete(name = "Transcription") +
  theme_classic() 
dev.off()
theme(legend.position = "bottom")
p <- plot_grid(m_to_d,b_to_d, ncol = 2)
ggsave("~/Desktop/supplement_final/log2_m_to_d_and_b_to_d.pdf", p, width = 12, height = 10)

############### Interaction tissue:treatment, how each tissue responds
## differently
pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig2_day_tiss_treat.pdf", width = 6, height = 10)
ggplot(deseq_selected_contigs_rev, aes(x = genes, y = day_tissue_treat_log2FoldChange, shape = day_tissue_treat_padj < 0.05)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = day_tissue_treat_log2FoldChange > 0 )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_point(aes(y = 3),
             position = position_dodge(0.9), 
             show.legend = FALSE) +
  scale_shape_manual(values = c(NA, 8)) +
  coord_flip()
dev.off()

pdf("~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/fig2_night_tiss_treat.pdf", width = 6, height = 10)
ggplot(deseq_selected_contigs_rev, aes(x = genes, y = night_tissue_treat_log2FoldChange, shape = night_tissue_treat_padj < 0.05)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = night_tissue_treat_log2FoldChange > 0 )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_point(aes(y = 3),
             position = position_dodge(0.9), 
             show.legend = FALSE) +
  scale_shape_manual(values = c(NA, 8)) +
  coord_flip()
dev.off()

####################### CCM GENES DE differently IN M VS M IN WW or DD (figure A)
### which genes are are high in differentially expressed across tissues only in WW or In D, but not in both
## of the list

ccm_in_d_only<-CCM_BS_vs_M_in_D[-which(CCM_BS_vs_M_in_D$V1 %in% CCM_BS_vs_M_in_WW$V1),]

ccm_in_ww_only<-CCM_BS_vs_M_in_WW[-which(CCM_BS_vs_M_in_WW$V1 %in% CCM_BS_vs_M_in_D$V1),]

# remove those with contigs already in the consisten M and B

#ccm_in_ww_d<-merge(ccm_in_ww[,c(2,3,5)],ccm_in_d[,c(2,3,5)],by='V1', all = TRUE)
ccm_in_d_only_1<-as.data.frame(ccm_in_d_only %>% group_by(V2) %>% top_n(1, log2FoldChange))
ccm_in_d_only_1<-ccm_in_d_only_1[-which(ccm_in_d_only_1$V2 %in% ccm_in_ww_d_1$V2.x),]

pdf("CCM_DE_M_vs_BS_only_in_D.pdf", width = 3.67, height = 2.20)
ggplot(ccm_in_d_only_1, aes(x = reorder(V2, -log2FoldChange), y = log2FoldChange)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = log2FoldChange > 0 )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()
dev.off()

## NOW WW
ccm_in_ww_only_1<-as.data.frame(ccm_in_ww_only %>% group_by(V2) %>% top_n(1, log2FoldChange))
ccm_in_ww_only_1<-ccm_in_ww_only_1[-which(ccm_in_ww_only_1$V2 %in% ccm_in_ww_d_1$V2.x),]

pdf("CCM_DE_M_vs_BS_only_in_WW.pdf", width = 3.67, height = 2.20)
ggplot(ccm_in_ww_only_1, aes(x = reorder(V2, -log2FoldChange), y = log2FoldChange)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = log2FoldChange > 0 )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()
dev.off()

######
#### NOW PLOT COMMON GENES DE IN BS AND M IN RESPONSE TO DROUGHT.
##############
### genes DE afte Drought in both M and BS
ccm_m_in_d<-CCM_DE_M_in_D[which(CCM_DE_M_in_D$V1 %in% CCM_DE_BS_in_D$V1),]

ccm_b_in_d<-CCM_DE_BS_in_D[which(CCM_DE_BS_in_D$V1 %in% CCM_DE_M_in_D$V1),]

ccm_b_and_m_in_d<-merge(ccm_b_in_d[,c(2,3,5)],ccm_m_in_d[,c(2,3,5)],by='V1', all = TRUE)

###
mtabTPM = melt(ccm_b_and_m_in_d[,-c(1,4)], id.vars="V2.x")
reorder(V2, -log2FoldChange)
# if several contigs, select the one with more change: 
library(dplyr)
ccm_b_and_m_in_d_1<-as.data.frame(ccm_b_and_m_in_d %>% group_by(V2.x) %>% top_n(1, log2FoldChange.x))

####

#####
data<-ccm_b_and_m_in_d_1[,-c(1,4)]
ggplot(data,
       aes(y = reorder(V2.x, -log2FoldChange.x))) +
  labs(x = "log2FoldChange_M_vs_BS",y= 'gene') +
  geom_segment(aes(x = data$"log2FoldChange.x",
                   y = reorder(V2.x, -log2FoldChange.x),
                   xend = data$"log2FoldChange.y",
                   yend = reorder(V2.x, -log2FoldChange.x)),
               size = 1) +
  geom_point(aes(x = data$"log2FoldChange.x", color = log2FoldChange.x > 0), size = 4, shape = 16) +
  geom_point(aes(x = data$"log2FoldChange.y", color = log2FoldChange.y > 0), size = 4, shape = 17) +
  scale_color_discrete(name = "tissue") +
  theme_classic() 
theme(legend.position = "bottom")
######################
###############
#######


####################### CCM GENES DE UNIQUELY IN M OR IN B in response to D

ccm_m_in_d_only<-CCM_DE_M_in_D[-which(CCM_DE_M_in_D$V1 %in% CCM_DE_BS_in_D$V1),]

ccm_b_in_d_only<-CCM_DE_BS_in_D[-which(CCM_DE_BS_in_D$V1 %in% CCM_DE_M_in_D$V1),]


#ccm_in_ww_d<-merge(ccm_in_ww[,c(2,3,5)],ccm_in_d[,c(2,3,5)],by='V1', all = TRUE)
ccm_m_in_d_only_1<-as.data.frame(ccm_m_in_d_only %>% group_by(V2) %>% top_n(1, log2FoldChange))
### if the contigs already are in the table of genes DE in both tissues, remove
ccm_m_in_d_only_1<-ccm_m_in_d_only_1[-which(ccm_m_in_d_only_1$V2 %in% ccm_b_and_m_in_d_1$V2.x),]

pdf("CCM_DE_in_B_only_in_D.pdf", width = 3.67, height = 2.20)
ggplot(ccm_b_in_d_only_1[-c(3,)], aes(x = reorder(V2, -log2FoldChange), y = log2FoldChange)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = log2FoldChange > 0 )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()
dev.off()

## NOW BS
#ccm_in_ww_d<-merge(ccm_in_ww[,c(2,3,5)],ccm_in_d[,c(2,3,5)],by='V1', all = TRUE)
ccm_b_in_d_only_1<-as.data.frame(ccm_b_in_d_only %>% group_by(V2) %>% top_n(1, log2FoldChange))
### if the contigs already are in the table of genes DE in both tissues, remove
ccm_b_in_d_only_1<-ccm_b_in_d_only_1[-which(ccm_b_in_d_only_1$V2 %in% ccm_b_and_m_in_d_1$V2.x),]

pdf("CCM_DE_in_B_only_in_D.pdf", width = 3.67, height = 2.20)
ggplot(ccm_b_in_d_only_1, aes(x = reorder(V2, -log2FoldChange), y = log2FoldChange)) + 
  theme_classic() + geom_bar(stat = "identity", aes(fill = log2FoldChange > 0 )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()
dev.off()

##########



ggplot(mtabTPM,aes(V2.x, value)) +
  geom_point(aes(color=year)) +
  geom_line(aes(group = paired))
#sig_dif_response_in_tissue<-as.data.frame(sig_dif_response_in_tissue)
#ccm_in_ww_d<-merge(ccm_in_ww_d,sig_dif_response_in_tissue, by =0, all=F)
to_plot<-(ccm_in_ww_d[,c(3,5)])
to_plot<-to_plot[order(to_plot$log2FoldChange),] 
# renae duplicated gene identities
to_plot$V2<-make.names(to_plot$V2,unique=T)

####



### how to know if they are of the same direction of change and 

###

##### more understanding: 
# https://support.bioconductor.org/p/92932/
#ds <- DESeqDataSetFromMatrix(countData=day_counts, colData=samples, design=~tisse + treatment + tisse:treatment)
