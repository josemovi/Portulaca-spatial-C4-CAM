
## FROM EXTRACT_DE_ccm_anotatated_genes_from_kallisto_results_3june2021.r

########
ct<-read.delim('~/Desktop/supplement_final/lines_plots_DEseq_table_ccm_manual_reduced.csv', sep = '\t')
ct2<-read.delim('~/Desktop/supplement_final/lines_plots_DEseq_table_ccm_manual_reduced.csv', sep = '\t')

#########
#######
row.names(ct)<-ct$indexes
########

tpm<-as.data.frame(so_norm_filt_matrix[row.names(so_norm_filt_matrix) %in% ct$unigenes,])
tpm$unigenes<-rownames(tpm)
ct$enzyme<-row.names(ct)
#### reduce counts table to list of selected genes
ct<-merge(ct,tpm,by = 'unigenes')
row.names(ct)<-ct$enzyme
## GET ONLY THE TPME DATA
exp_tab<-as.data.frame(t(ct[,c(16:43)]))
#### make columins with plant id, group of samples, time, tissue
g1<-(str_split_fixed(as.list(row.names(exp_tab)), "_",2)[,2])
g2<-str_split_fixed(g1, "-",2)[,2]
g2<-gsub("DD", "D", g2)

exp_tab$group<-g2
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
library(reshape2)
mtabTPM = melt(exp_tab, id.vars=c("group","treatment","time","tissue","timetissue","treatmenttissue"))

######### WITH IQR: 
df_iqr<-as.data.frame(unclass(t(medIQR(value ~ group : variable , data = mtabTPM))))

#### make columins with plant id, group of samples, time, tissue
df_iqr$group<-(str_split_fixed(as.list(row.names(df_iqr)), "[.]",2)[,1])
df_iqr$enzyme<-(str_split_fixed(as.list(row.names(df_iqr)), "[.]",2)[,2])

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
#################### PLOT BY MODULE: 
df_iqr_new <- df_iqr                              # Replicate data
df_iqr_new$enzyme <- factor(df_iqr_new$enzyme,      # Reordering group factor levels
                         levels = ct2$indexes)

#enz_plots = list()
#for (val in 1:length(levels(ct$Pathway))) {
#  enz<-rownames(ct[grep(levels(ct$Pathway)[val],ct$Pathway),])
#  reduced_df<-df_iqr[which(df_iqr$enzyme %in% enz),]
#  enz_plots[[val]]<-ggplot(reduced_df, aes(x=time, y=Median, group=treatmenttissue, color=tissue)) + 
p<-ggplot(df_iqr_new, aes(x=time, y=Median, group=treatmenttissue, color=tissue)) + 
  geom_line(aes(linetype=treatment, color=tissue), size = 1) +
  geom_point()+
#  ggtitle(levels(ct$Pathway)[val]) +
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.1) +
  theme_classic() +
  scale_color_manual(values=cols) +
  scale_linetype_manual(values=lin) +
  theme(strip.text = element_text(size=10)) +
  facet_wrap( ~ enzyme, scales="free", ncol = 5)
#}
#dev.off()
#### PLOT
#nrow(ct[grep(levels(ct$Pathway)[1],ct$Pathway),])

#pdf("~/Desktop/supplement_final/lines_plots/Carboxylation.pdf", height = 4, width = 7)
#enz_plots[[2]]
#dev.off()
#ggsave(
#  filename = "~/Desktop/supplement_final/line_plots.pdf", 
#  plot = marrangeGrob(enz_plots, nrow=1, ncol=1),  
#  width = 4, height = 9)

#p <- plot_grid(enz_plots[[1]],enz_plots[[2]], enz_plots[[3]],enz_plots[[4]],enz_plots[[5]], ncol = 1, nrow = 5)
#grid.arrange(enz_plots[[1]],enz_plots[[2]],nrow = 2)
pdf('~/Desktop/supplement_final/line_plots.pdf', height = 35, width = 10)
p
dev.off()

###################################################
#####################################################
#####################################

### plot heatmap
head(df_iqr)
head(df_iqr[,c(1,4,5)])
df_iqr_dcast<-dcast(df_iqr[,c(1,4,5)], enzyme ~ group, value.var = "Median")
df_iqr_dcast$module<-(str_split_fixed(as.list(row.names(df_iqr_dcast)), " ",2)[,1])
out[order(as.numeric(as.character(out$V3))), ]
df_iqr_dcast<-df_iqr_dcast[order(as.numeric(as.character(df_iqr_dcast$module))), ]
#####
row.names(df_iqr_dcast)<-df_iqr_dcast[,1]
db<-as.matrix(df_iqr_dcast[,-c(1,10)])
####
names(df_iqr_dcast)
library(pheatmap)
db_relative<-round(db/rowSums(db), 2)*100
db2<-db
head(db)
rownames(db)<-paste(rownames(db),round(apply(db, 1, max), digits = 1),sep = ' ')
db_relative_to_max <- t(apply(db, 1, function(x) x/max(x) * 100))
head(db_relative_to_max)
###
#pheatmap(db_relative[c(1:102),c(1,5,2,6,3,7,4,8)],
pdf('../kallisto_mRNA_sleuth/sleuth_gene_mode_13may21/heatmap_sleuth_significant_DE_representative_contig_per_expression_class.pdf', height = 15, width = 6)
pheatmap(db_relative_to_max[c(1:102),c(6,2,5,1,8,4,7,3)],
         cluster_cols = F,cluster_rows = T,
         colorRampPalette(c("white", "lightblue","lightpink","firebrick3"))(4),
#         gaps_row = c(10, 14), 
         gaps_col = 4)
dev.off()
#pheatmap(db_relative[c(1:102),c(1,5,3,7,2,6,4,8)],cluster_cols = F,cluster_rows = F)
###