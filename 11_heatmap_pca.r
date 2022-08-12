so_norm_filt_matrix
### 
library(sleuth)
log_sleuth_tpm<-log_transform(so_norm_filt_matrix, offset = 0.5)
var_genes <- apply(log_sleuth_tpm, 1, var)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- log_sleuth_tpm[select_var,]
dim(highly_variable_lcpm)
#########
library(pheatmap)
pdf("~/Desktop/supplement_final/libraries_heatmap.pdf", width = 6, height = 4)
db<-as.matrix(highly_variable_lcpm)
pheatmap(db, show_rownames = F)
dev.off()
##################################

##################################################################
###################################################################
### pca WITH CIRCLES ########################################
library(factoextra)
pca_tpm<-t(highly_variable_lcpm)
pca_tpm<-(log_sleuth_tpm)
pca_tpm[pca_tpm == 0] <- 0.1
res.pca <- prcomp(pca_tpm, scale = TRUE)
### extract results from variables: 
var <- get_pca_ind(res.pca)
head(var$contrib)
#####
p <- pca(log_sleuth_tpm, removeVar = 0.1)

library(factoextra)
head(pca_tpm)
groups <- as.factor(g2)
cols <- c('WW.7.M'='green','WW.7.BS'='darkgreen',
          'WW.19.M'='violet','WW.19.BS'='darkviolet',          
          'WW.23.M'='blue','WW.23.BS'='darkblue',
          'WW.3.M'='grey','WW.3.BS'='black',
          'D.7.M'='yellow','D.7.BS'='khaki',
          'D.15.M'='coral','D.15.M'='coral4',
          'D.23.BS'='brown','D.23.M'='darkorange',
          'D.3.M'='red','D.3.M'='red4')
col_list<-as.data.frame(col.sample)[,1]
pdf("~/posdoc/florida_trans_spatial_analyses-dir/figures_spatial_transcriptomics_26March21/FigSX.log_PCA.pdf", width = 6, height = 4)
fviz_pca_biplot(res.pca,
#             col.ind = groups, # color by groups
             #             palette = c('khaki',"black","darkgreen","darkblue","orange","brown","lightgreen","lightblue","yellow","cyan","red",
             #                         "violet",'grey',"green",'coral','blue'),
#             palette = col_list,
#             addEllipses = TRUE, # Concentration ellipses
#             ellipse.type = "confidence",
#             legend.title = "Groups",
             repel = TRUE
)
dev.off()
