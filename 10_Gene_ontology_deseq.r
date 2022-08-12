
#############################
#############################
### GENE ONTOLOGY
#############################
BiocManager::install("topGO")
library(topGO)
library(ALL)
# in bash farnam: 
#https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-GOSeq
#Extract GO assignments per gene
cd /home/jjm247/project/florida_transcriptome_dir/Assembled-transcriptomes
~/programs/Trinotate-Trinotate-v3.2.1/util/extract_GO_assignments_from_Trinotate_xls.pl \
--Trinotate_xls trinotate_annotation_report.xls \
-G --include_ancestral_terms \
> florida_go_annotations.txt

################ copy to mac: 
go_annot<-read.delim('~/posdoc/florida_trans_spatial_analyses-dir/florida_go_annotations.txt', header = FALSE)
names(go_annot)<-c('target_id','GO_groups')


#http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
#####Reading in GO annotations for the genes
geneID2GO <- readMappings(file = "~/posdoc/florida_trans_spatial_analyses-dir/florida_go_annotations.txt")  
#####
geneUniverse <- names(geneID2GO) 
## genes of interest
#treatmentD_sig$unigenes,
#Tissue = tisseM_sig$unigenes,
#Time = time7_sig$unigenes,
#Individual = plantPO2_sig$unigenes
genesOfInterest <- plantPO2_sig$unigenes
#genesOfInterest <- time7_sig$unigenes
#genesOfInterest <- tisseM_sig$unigenes
#genesOfInterest <- treatmentD_sig$unigene
###### Then we need to tell TopGO where these interesting genes appear in the 'geneUniverse' vector:
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
#####

names(geneList) <- geneUniverse
####Putting the data together into an R object
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
#myGOdata 
##########################
resultFisher <- runTest(myGOdata, algorithm="elim", statistic="fisher") 

allRes_Plant <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 50)
#allRes_Time <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 50)
#allRes_Tissue <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 50)
#allRes_Treatment <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 50)

##########################
###########PLOT NUMBER OF GENES ########
##############################
pdf("Fig2_B.GO_plant.pdf", width = 6, height = 8)
p<-ggplot(data=allRes_Time, aes(x=reorder(paste0(allRes_Time$Term,allRes_Time$GO.ID,sep ='-'), Significant), y= Significant)) +
  #ggplot(data=allRes, aes(x=reorder(Term, -Significant), y= Significant)) +
  labs(x = "GO Term", y = "Unigene counts") +
  geom_bar(stat="identity",  fill="steelblue")
p + coord_flip() +
  theme_minimal()
dev.off()
#########################
#############################
##########################
###########PLOT ENRICHMENT OF UP AND DOWN GENES ########
##############################
nrow(resultFisher)
### make plot similar to this one: https://www.researchgate.net/figure/Gene-Ontology-GO-Enrichment-Forest-Plot-of-postischemia-compared-to-baseline-in_fig4_270910808
######################### 
allRes<-allRes_Treatment
variable_qval<-'DvsWW_padj'
variable_b<-'DvsWW_log2FoldChange'
list_of_prop_expr <- list()
#temp<-data.frame('-1'=0,'1'=0)
#names(temp)<-c('-1','1')
for (i in 1:length(allRes$GO.ID)){
  test<-go_annot[grep(allRes$GO.ID[i],go_annot$GO_groups),]
  #  test<-go_annot[grep('GO:0007035',go_annot$GO_groups),]
  # significant
  test <-deseq_all_factors_tab[deseq_all_factors_tab$unigenes %in% test$target_id, ]
  test2 <- test[test[,variable_qval] < 0.05,]
  ## keep only one contig per unigene 
  test2<-test2[!duplicated(test2$unigenes),]
  ####
  test2<-sign(test2[,variable_b])
  ### check how many of each sign
  a<-sum(test2[test2==-1], na.rm=TRUE)*-1
  b<-sum(test2[test2==1], na.rm=TRUE)
  tb1<-prop.table(c(a,b))
  list_of_prop_expr[[i]]<-append(tb1,allRes$GO.ID[i],0)
}
df <- as.data.frame(do.call("rbind",list_of_prop_expr))
#names(df)<-c('GO.ID','BS','M')
names(df)<-c('GO.ID','WW','D')
#names(df)<-c('GO.ID','Night','Day')

df<-merge(df,allRes, by ='GO.ID')
#df$Night<-as.numeric(as.character(df$Night))
df$D<-as.numeric(as.character(df$D))
df$WW<-as.numeric(as.character(df$WW)) * -1
####
df2<-df[with(df, order(WW)),]
row.names(df2)<-c(1:nrow(df2))
df2<-df2[-c(1:7),]
row.names(df2)<-c(1:nrow(df2))
#df3<-df2[with(df2, order(-BS)),]
## melt
#library(reshape2)
mdf<-melt(df2[,c(1:4)], id.vars=c('Term','GO.ID'))
names(mdf)[c(3,4)]<-c('Treatment','Percent_genes')
#####################################################
pdf("~/posdoc/florida_trans_spatial_analyses-dir/supplement_final/figures/Fig_S.2B_GO_treatment.pdf", width = 6, height = 7)
ggplot(data=mdf, aes(x=reorder(paste0(mdf$Term,mdf$GO.ID,sep ='-'),Percent_genes), y= Percent_genes, fill = Treatment)) +
  #ggplot(data=allRes, aes(x=reorder(Term, -Significant), y= Significant)) +
  labs(x = "GO Term", y = "Percentage unigenes") +
  geom_bar(stat="identity",  position="identity") +
  theme(text = element_text(size=0.5)) +
  coord_flip() +
  theme_minimal()
dev.off()
###########################
##########################

