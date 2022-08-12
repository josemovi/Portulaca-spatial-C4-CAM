### extract intersect of DE genes: 
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


x <- list(
  Treatment = treatmentD_sig$unigenes,
  Tissue = tisseM_sig$unigenes,
  Time = time7_sig$unigenes,
  Individual = plantPO2_sig$unigenes
)

# get intersection of the 3 contrasts
DE_all_groups<-Reduce(intersect, x[-4])
trinotate<-read.delim('~/posdoc/florida_trans_spatial_analyses-dir/trinotate_annotation_report.xls')
trinotate$unigenes = sub("_i.*", "", as.character(trinotate$transcript_id))
head(DE_all_groups)
DE_all_groups<-trinotate[trinotate$unigenes %in% DE_all_groups,][c('unigenes','sprot_Top_BLASTP_hit')]
DE_all_groups_ccm<-comb_annot[comb_annot$unigenes %in% DE_all_groups$unigenes,]
DE_m_t_T<-DE_all_groups[!(DE_all_groups$unigenes %in% DE_all_groups_ccm$unigenes),][c('unigenes','sprot_Top_BLASTP_hit')]

DE_m_t_T<-DE_m_t_T[!duplicated(DE_m_t_T$unigenes),]
write.table(DE_m_t_T,'~/Desktop/supplement_final/DE_genes_annot_intersect_all_factors1.csv')

#########
#genes<-rbind(genesM,genesD,genestime)
############# ADD INTERSECTION ALL OF THEM
#genes_intersect<-DE_m_t_T[grep('Scarecrow|MPI1|Potassium transporter|ATPase|PUMP|DTC|ACOX|LETM|AAPC|NDH|PLASMODESMATA|UDP|P2A01|ERD',DE_m_t_T$sprot_Top_BLASTP_hit),]
#genes_intersect$Pathway<-'intersect'
#head(genes_intersect)
#genes<-rbind(genes,genes_intersect)

#DE_m_t_T$enzyme = sub("_.*", "", as.character(DE_m_t_T$sprot_Top_BLASTP_hit))
#DE_m_t_T<-DE_m_t_T[duplicated(DE_m_t_T$enzyme),]
#DE_m_t_T[,c(1,3)]
#################PLOTTT
tp2<-as.data.frame(so_norm_filt_matrix[row.names(so_norm_filt_matrix) %in% DE_m_t_T$unigenes,])

tp2$unigenes<-rownames(tp2)
reduced_genes<-read.csv('~/Desktop/supplement_final/DE_genes_annot_intersect_all_factors_manually_reduces.csv',sep = ' ')
tp2<-merge(DE_m_t_T,tp2,by = 'unigenes')
tp2<-merge(reduced_genes,tp2,by = 'unigenes')
tp2$enzyme = sub("_.*", "", as.character(tp2$sprot_Top_BLASTP_hit))
row.names(tp2)<-paste(tp2$unigenes,tp2$enzyme,sep="_")

## GET ONLY THE TPME DATA
exp_tab<-as.data.frame(t(tp2[,c(4:30)]))
exp_tab<-round(exp_tab,digits = 2)
#### make columins with plant id, group of samples, time, tissue
g1<-(str_split_fixed(as.list(row.names(exp_tab)), "_",2)[,2])
g2<-str_split_fixed(g1, "-",2)[,2]
g2<-gsub("DD", "D", g2)

exp_tab$group<-g2
## plant
#pca_tpm_startch$sample<-row.names(pca_tpm_startch)
#pca_tpm_startch$plant<-str_split_fixed(g1, "[.]",2)[,1]
# treatment
tr<-(str_split_fixed(as.list(exp_tab$group), "-",3)[,1])
exp_tab$treatment<-tr
time<-(str_split_fixed(as.list(exp_tab$group), "-",3)[,2])
exp_tab$time<-time
## tissue
tissue<-(str_split_fixed(as.list(exp_tab$group), "-",3)[,3])
exp_tab$tissue<-tissue
## time tissue
timetissue<-(str_split_fixed(as.list(exp_tab$group), "-",2)[,2])
exp_tab$timetissue<-timetissue
##tissue treatment
treatmenttissue<-as.data.frame(str_split_fixed(as.list(exp_tab$group), "-",3)[,c(1,3)])
exp_tab$treatmenttissue<-paste(treatmenttissue$V1,treatmenttissue$V2,sep="_")

head(exp_tab)

#library(data.table)
#install.packages("summarize", repos="http://R-Forge.R-project.org")
library(summarize)
#PO1<-exp_tab[grep('PO1',rownames(exp_tab)),]

mtabTPM = melt(exp_tab, id.vars=c("group","treatment","time","tissue","timetissue","treatmenttissue"))

######### WITH IQR: 
df_iqr<-as.data.frame(unclass(t(medIQR(value ~ group : variable , data = mtabTPM))))

#### make columins with plant id, group of samples, time, tissue
df_iqr$group<-(str_split_fixed(as.list(row.names(df_iqr)), "[.]",2)[,1])
df_iqr$enzyme<-(str_split_fixed(as.list(row.names(df_iqr)), "[.]",2)[,2])
df_iqr$gene<-(str_split_fixed(as.list(row.names(df_iqr)), "_",5)[,5])

tr<-(str_split_fixed(as.list(df_iqr$group), "-",3)[,1])
df_iqr$treatment<-tr
time<-(str_split_fixed(as.list(df_iqr$group), "-",3)[,2])
df_iqr$time<-time
## tissue
tissue<-(str_split_fixed(as.list(df_iqr$group), "-",3)[,3])
df_iqr$tissue<-tissue
## time tissue
timetissue<-(str_split_fixed(as.list(df_iqr$group), "-",2)[,2])
df_iqr$timetissue<-timetissue
##tissue treatment
treatmenttissue<-as.data.frame(str_split_fixed(as.list(df_iqr$group), "-",3)[,c(1,3)])
df_iqr$treatmenttissue<-paste(treatmenttissue$V1,treatmenttissue$V2,sep="_")


#### colors 
cols <- c('D_M'='red','D_BS'='red',
          'WW_BS'='black','WW_M'='black','D'='red','WW'='black','BS'='red','M'='black')
#cols <- c('D_M'='black','D_BS'='green',
#          'WW_BS'='green','WW_M'='black','D'='red','WW'='black')
lin <- c('D_M'='dashed','D_BS'='solid',
         'WW_BS'='solid','WW_M'='dashed','BS'='solid','M'='dashed','WW'='solid','D'='dashed')
####################
height = 35, width = 10)
pdf('~/Desktop/supplement_final/Fig_reduced_intersect_DE.pdf', height = 20, width = 10)
ggplot(df_iqr, aes(x=time, y=Median, group=treatmenttissue, color=tissue)) + 
  geom_line(aes(linetype=treatment, color=tissue), size = 1) +
  geom_point()+
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.1) +
  theme_classic() +
  scale_color_manual(values=cols) +
  scale_linetype_manual(values=lin) +
  theme(strip.text = element_text(size=10)) +
  facet_wrap( ~ enzyme, scales="free", ncol = 5) +
  theme(strip.text.x = element_text(size = 4))
dev.off()
######
