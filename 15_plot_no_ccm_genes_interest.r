## NO CCM GENES OF INTEREST
noccm<-read.csv('~/Desktop/supplement_final/list_no_ccm_genes.csv',header = F)

deseq_interst<-read.csv('~/Desktop/supplement_final/DE_genes_annot_intersect_all_factors1.csv', sep = ' ')
sleuth_interst<-read.csv('~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/sleuth_gene_mode_13may21/DE_genes_annot_intersect_all_factors_7june21.csv', sep = ' ')
#View(ds_noccm_trino)

ds_noccm_trino<-deseq_interst[grep(paste(noccm$V1,collapse="|"),deseq_interst$sprot_Top_BLASTP_hit),]
sl_noccm_trino<-sleuth_interst[grep(paste(noccm$V1,collapse="|"),sleuth_interst$sprot_Top_BLASTP_hit),]
others<-trinotate[grep("EXPA2",trinotate$sprot_Top_BLASTP_hit),]
noccm_trino<-rbind(ds_noccm_trino,sl_noccm_trino,others[,c('unigenes','sprot_Top_BLASTP_hit')])
noccm_trino<-noccm_trino[!duplicated(noccm_trino$unigenes),]
#########
### to report: 
noccm_trino$Description = sub(";.*", "", as.character(noccm_trino$sprot_Top_BLASTP_hit))
noccm_trino$Description = sub(".*=", "", as.character(noccm_trino$Description))
noccm_trino$Uniprot_gene = sub("_.*", "", as.character(noccm_trino$sprot_Top_BLASTP_hit))
write.csv(noccm_trino,'~/Desktop/supplement_final/table_for_selected_noccm.csv')
# get intersection of the 3 contrasts
#DE_all_groups<-Reduce(intersect, x[-4])
#trinotate<-read.delim('~/posdoc/florida_trans_spatial_analyses-dir/trinotate_annotation_report.xls')
#trinotate$unigenes = sub("_i.*", "", as.character(trinotate$transcript_id))

tp2<-as.data.frame(so_norm_filt_matrix[row.names(so_norm_filt_matrix) %in% noccm_trino$unigenes,])

tp2$unigenes<-rownames(tp2)

noccm_trino<-noccm_trino[!duplicated(noccm_trino$unigenes),]
tp2<-merge(noccm_trino,tp2,by = 'unigenes')
tp2$enzyme = sub("_.*", "", as.character(tp2$sprot_Top_BLASTP_hit))

row.names(tp2)<-paste(tp2$unigenes,tp2$enzyme,sep="_")

## GET ONLY THE TPME DATA
exp_tab<-as.data.frame(t(tp2[,c(3:30)]))
#exp_tab<-round(exp_tab,digits = 2)
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
## factor to sort genes alphabetically

df_iqr$gene <- factor(df_iqr$gene, levels = unique(sort(df_iqr$gene)))


#### colors 
cols <- c('D_M'='red','D_BS'='red',
          'WW_BS'='black','WW_M'='black','D'='red','WW'='black','BS'='red','M'='black')
#cols <- c('D_M'='black','D_BS'='green',
#          'WW_BS'='green','WW_M'='black','D'='red','WW'='black')
lin <- c('D_M'='dashed','D_BS'='solid',
         'WW_BS'='solid','WW_M'='dashed','BS'='solid','M'='dashed','WW'='solid','D'='dashed')
####################

pdf('~/Desktop/supplement_final/Fig_selected_no_ccm_DE.pdf', height = 9, width = 10)
ggplot(df_iqr, aes(x=time, y=Median, group=treatmenttissue, color=tissue)) + 
  geom_line(aes(linetype=treatment, color=tissue), size = 1) +
  geom_point()+
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.1) +
  theme_classic() +
  scale_color_manual(values=cols) +
  scale_linetype_manual(values=lin) +
  theme(strip.text = element_text(size=10)) +
  facet_wrap( ~ gene, scales="free", ncol = 5) +
  theme(strip.text.x = element_text(size = 6))
dev.off()
######