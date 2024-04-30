####
## PCA LASI-DAD 1000G
####

# load tidyverse package
library(tidyverse)

###
##Data
###
#filename
file_vec="LASI-DAD2620_1000G_EUR_EAS_SAS_pruned05_maf0.05_chr1-22_Fig1B.txt"
file_val="LASI-DAD2620_1000G_EUR_EAS_SAS_pruned05_maf0.05_chr1-22_Fig1B_eval.txt"
#load data
pca_1000G_LASI_DAD <- read_table(file_vec, col_names = T,na="NA")
eigenval <- scan(file_val)
pve <- data.frame(PC = 1:10, pve = eigenval[1:10]/sum(eigenval)*100)

###
##Plots in function of Region
###
#Eigenvalues
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

##Plot in function of Region of origin
b <- ggplot(filter(pca_1000G_LASI_DAD, Region!="SAS"), aes(PC1, PC2, col=Region,pch=Region)) + geom_point(size = 2) +
  scale_color_manual(breaks = c("North", "North-East","East","Central","West","South","Other"), 
                     labels = c("North", "North-East","East","Central","West","South","Other"),
                     values=c("#f58231","#4363d8","#42d4f4","#e6194B","#ffe119","#3cb44b","#a9a9a9",NA,NA,NA), name="LASI-DAD", na.value="lightgrey")+
  scale_shape_manual(breaks = c("EUR","EAS"), 
                     labels = c("European","East Asian"),
                     values=c("EUR"=18,"EAS"=17), name="1000G", na.value=20)+
  guides(colour = guide_legend(override.aes = list(size=7)), shape = guide_legend(override.aes = list(size=7, col="lightgrey")))
b <- b + coord_equal() + theme_bw() +
  xlab(sprintf("PC1 (%.2f%s)",pve$pve[1],"%"))+ylab(sprintf("PC2 (%.2f%s)",pve$pve[2],"%"))+
  theme_light(base_size = 20)
b
