### create supplemental table for main text

ct<-read.delim('~/Desktop/supplement_final/lines_plots_DEseq_table_ccm_manual_reduced.csv', sep = '\t')
deseq_all_factors_tab_trinot_2<-read.csv('~/posdoc/florida_trans_spatial_analyses-dir/kallisto_mRNA_sleuth/DEseq_august21/table_S4_trinotate_deseq_v1',sep = "\t")
######### merge means
head(selected_contigs)

#####

table_unigenes_sup<-merge(deseq_all_factors_tab_trinot_2,ct, by= 'unigenes')
table_unigenes_sup<-table_unigenes_sup[match(ct$unigenes, table_unigenes_sup$unigenes),]

nrow(table_unigenes_sup)
write.csv(table_unigenes_sup,'~/Desktop/supplement_final/supplement_table_genes.csv')
