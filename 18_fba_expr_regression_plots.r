fba<-read.csv('~/Desktop/figures_lcm_final/EXPFBA_gene_corrected.csv')
fba$Enzyme <- factor(fba$Enzyme, levels = unique(fba$Enzyme))

cols5 <- c('BS'='blue4','M'='cyan')
fba$Tissue<-gsub('.*_','',fba$Location)
my.formula <- EXP ~ FBA
set.seed(42)
library(ggpubr)
require("ggrepel")
ggplot(fba,aes(x = FBA, y = EXP)) +
  geom_point(size = 3, aes(shape=Treatment, colour = Tissue)) + 
#  geom_text(label = fba$Time) +
   geom_text_repel(label = fba$Time) +
#  geom_label_repel(aes(label = fba$Time, 
#                      fill = factor(fba$Time)), color = 'white') +
  scale_colour_manual(values = cols5) +
  scale_shape_manual(values=c(8, 17)) +
  geom_smooth(method='lm') +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.y= 1.4, label.x = 0.9, vjust = 1, hjust = 1.1, size = 2.5) +
  facet_wrap( ~ Enzyme, scales="free", ncol = 4) +
  scale_y_continuous(limits = c(-0.2,1.4)) +
  scale_x_continuous(limits = c(-0.2,1.4)) +
  theme_bw()

