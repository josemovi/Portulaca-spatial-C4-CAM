### R commands to to test for dial acid accumulation and plot results. 

### Author: J.J. Moreno-Villena. Supplement Code for 'Spatial resolution of an integrated C4+CAM photosynthetic metabolism'

## PLOT RESULTS
# Read table
titra<-read.table('titrations_data.csv', sep = '\t', header = T)

# Summary function from: ttp://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# Sumary of the table
df2 <- data_summary(titra, varname="uequiH", 
                      groupnames=c("group")) #"clone"))
                     
# Split group label into treatment and time labels 
df2$treatment<-(str_split_fixed(df2$group, "-",2)[,1])
df2$time<-(str_split_fixed(df2$group, "-",2)[,2])

#Factorize groups
df2$time = factor(df2$time,levels=c("19:00", "07:00"))
df2$treatment = factor(df2$treatment,levels=c("Watered", "Drought"))

# Plot barpltos with ggplot
library(ggplot2)
p<-ggplot(df2, aes(x=time, y=uequiH, fill=treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=uequiH-sd, ymax=uequiH+sd), width=.2,
                position=position_dodge(.9)) +
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))
# Save 
pdf('titrations.pdf', width=5,height=4)
p
dev.off()

## T-test for diel acidity accumulation in well watered samples
t.test(titra[grep('Watered-07:00',titra$group),14],titra[grep('Watered-19:00',titra$group),14], alternative = "greater", var.equal = FALSE)

## T-test for diel acidity accumulation in droughted sample
t.test(titra[grep('Drought-07:00',titra$group),14],titra[grep('Drought-19:00',titra$group),14], alternative = "greater", var.equal = FALSE)
