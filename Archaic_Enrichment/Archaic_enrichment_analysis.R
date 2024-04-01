####
## PBS analysis
####
library(ggplot2)
library(ggpattern)
library(dplyr)
library(tidyr)

setwd("/LASI-DAD/Archaic_Enrichment/")

##import data
data_ND_LASIDAD_EUR_EAS = read.table('archaic_freq_perpop_LASIDAD_1000G_EUROPE_EASTASIA_0.8_PHASED.txt', header = T, sep = '\t')


###chromosome length in hg38
chromosomes = data.frame(chrom = c(1:22),
                         starts = rep(0,22),
                         ends = c(248957*1000,242194*1000,198296*1000,190215*1000,181539*1000,170806*1000,159346*1000,145139*1000,138395*1000,133798*1000,135087*1000,133276*1000,114365*1000,107044*1000,101992*1000,90339*1000,83258*1000,80374*1000,58618*1000,64445*1000,46710*1000,50819*1000)) %>%
  mutate(chrom = factor(chrom, c(1:22)) ) %>%
  mutate(labels = paste0('chr ', chrom, ' '))

##import gaps
gaps = read.table('callability.bed', header = F,  col.names = c('chrom','start','end')) %>%
  filter(end - start > 50000, chrom != 'X') %>%
  mutate(chrom = factor(chrom, c(1:22)) )



####
##Distribution frequency Archaic in LASIDAD
####
average_NEA_LASIDAD=mean(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_NEA_freq)
sd_NEA_LASIDAD=sd(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_NEA_freq)
lim_enriched_NEA=average_NEA_LASIDAD+2*sd_NEA_LASIDAD
#top 5% callable
percentile_95N=quantile(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_NEA_freq,0.95)
percentile_99N=quantile(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_NEA_freq,0.99)

##denisovan
average_DEN_LASIDAD=mean(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_DEN_freq)
sd_DEN_LASIDAD=sd(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_DEN_freq)
lim_enriched_DEN=average_DEN_LASIDAD+2*sd_DEN_LASIDAD
#top 5% callable
percentile_95D=quantile(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_DEN_freq,0.95)
percentile_99D=quantile(filter(data_ND_LASIDAD_EUR_EAS,called>0.5)$LASIDAD_DEN_freq,0.99)



#####
## Enrichment in LASI-DAD
#####
##group and window
window_size=4e4

data_mutate_99freq = data_ND_LASIDAD_EUR_EAS %>%
  mutate(rounded_start = start - start%%window_size) %>%
  group_by(chrom, rounded_start) %>%
  summarize(
    mean_NEA_freq = mean(LASIDAD_NEA_freq), 
    mean_DEN_freq = mean(LASIDAD_DEN_freq),
    mean_called = mean(called) )%>%
  ungroup() %>%
  gather(archaic, freq, mean_NEA_freq:mean_DEN_freq) %>%
  mutate(Yaxis = ifelse(archaic == 'mean_NEA_freq', 1,2 ), 
         top99 = ifelse((archaic == 'mean_NEA_freq' & freq > lim_enriched_NEA & mean_called>0.5) | (archaic == 'mean_DEN_freq' & freq > lim_enriched_DEN & mean_called>0.5),"1","0" )) 


##add desert
chrom_desert_N=c(3,5,7,8,8,18)
start_desert_N=c(76800000,83000000,106700000,52400000,107800000,32300000)
end_desert_N=c(90400000,94300000,128200000,65400000,121500000,46300000)
Yaxis=rep(1,length(chrom_desert_N))
desert_N<- data.frame(chrom_desert_N,start_desert_N,end_desert_N,Yaxis)
colnames(desert_N)[1]="chrom"
colnames(desert_N)=c("chrom","start_desert","end_desert","Yaxis")

chrom_desert_D=c(1,2,2,3,3,3,5,7,10,12,13,14,14)
start_desert_D=c(108400000,98400000,147700000,77000000,93700000,130600000,58800000,81700000,96400000,78900000,52300000,53400000,76300000)
end_desert_D=c(118800000,108700000,158000000,90400000,104600000, 143600000,69500000,97900000,106800000,93800000,
               63000000,64500000,90500000)
Yaxis=rep(2,length(chrom_desert_D))
desert_D <- data.frame(chrom_desert_D,start_desert_D,end_desert_D,Yaxis)
colnames(desert_D)[1]="chrom"
colnames(desert_D)=c("chrom","start_desert","end_desert","Yaxis")

archaic_deserts<- rbind(desert_N,desert_D) 

##to plot Neanderthal and Denisovan
Individuals = c(1,2)


###Plot 99%
data_mutate_99freq %>%
  mutate(chrom = factor(chrom, c(1:22) )) %>%
  
  ggplot() +
  
  facet_grid(chrom~., switch = 'y') +
  
  # chromosomes positions
  geom_rect(data = chromosomes, aes(xmin = starts/1e6, xmax = ends/1e6, ymin = 0, ymax = length(Individuals) + 0.01), fill = 'white', color = 'white') +
  geom_text(data = chromosomes, aes(x = 0, y = (length(Individuals))/3, label = labels), hjust = 1) +
  geom_rect(data = gaps, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.02, ymax = length(Individuals) + 0.005), fill = 'white') +
  geom_rect(data = chromosomes, aes(xmin = starts/1e6, xmax = ends/1e6, ymin = 0, ymax = length(Individuals) + 0.005), fill = NA, color = 'white') +
  
  # Frequency in windows
  geom_rect(aes(xmin = rounded_start/1e6, xmax = (rounded_start + window_size)/1e6, ymin = Yaxis-0.99, ymax = Yaxis-0.01, alpha=top99, fill=as.factor(Yaxis))) +
  
  ##add gaps
  geom_rect(data = gaps, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.02, ymax = length(Individuals) + 0.005), fill = 'white') +
  
  ##add desert N
  geom_rect_pattern(data = archaic_deserts, aes(xmin = start_desert/1e6, xmax = end_desert/1e6, ymin = Yaxis-0.99, ymax = Yaxis + 0.005 , pattern="Desert",fill = as.factor(Yaxis)),
                    col = "black",
                    pattern_spacing = 0.2, pattern_density = 0.2,
                    pattern_key_scale_factor = 0.2) +
  
  scale_fill_manual(values = c("#56B4E9","#009E73"),breaks = c("1","2"), labels=c("Neanderthal","Denisovan"),name="Ancestry") +
  scale_alpha_manual(breaks=c("0","1"), labels=c("No","Yes"), values=c(0.1,0.9), name="Enriched (>μ+2σ)")+
  scale_pattern_manual(values = c("Desert" = "stripe"), name="")+
  theme_bw(base_size=15) +
  labs(x = 'genomic position (Mb)') +
  theme(strip.background = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_blank(), #element_text(angle = 180)
        legend.position = c(0.8, 0.2),
        panel.spacing.y=unit(0, "lines"))+
  guides(fill = guide_legend(override.aes = list(pattern="none")), pattern = guide_legend(override.aes = list(fill="white", pattern="stripe")))




#####
## PBS in LASI-DAD compared to 1000G EUR and EAS
#####

##
# Functions definition for Fst
##


##is nan for dataframe
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))


Fst_2freq <- function(dataset,col_a,col_b){
  p_est=((dataset[col_a]/100+dataset[col_b]/100)/2)
  q_est=(((1-dataset[col_a]/100)+(1-dataset[col_b]/100))/2)
  pq_obs=(((dataset[col_a]/100)*(1-dataset[col_a]/100))+((dataset[col_b]/100)*(1-dataset[col_b]/100)))/2
  n=p_est*q_est-pq_obs
  d=p_est*q_est
  Fst=n/d
  Fst[is.nan(Fst)] <- 0
  return(Fst)
}

##Get column number from dataframe
col_NEA_EAS=which( colnames(data_ND_LASIDAD_EUR_EAS)=="EAS_NEA_freq" )
col_DEN_EAS=which( colnames(data_ND_LASIDAD_EUR_EAS)=="EAS_DEN_freq" )
col_NEA_EUR=which( colnames(data_ND_LASIDAD_EUR_EAS)=="EUR_NEA_freq" )
col_DEN_EUR=which( colnames(data_ND_LASIDAD_EUR_EAS)=="EUR_DEN_freq" )
col_NEA_LASIDAD=which( colnames(data_ND_LASIDAD_EUR_EAS)=="LASIDAD_NEA_freq" )
col_DEN_LASIDAD=which( colnames(data_ND_LASIDAD_EUR_EAS)=="LASIDAD_DEN_freq" )

###Neanderthal
##Fst
data_ND_LASIDAD_EUR_EAS$Fst_NEA_EAS_EUR<-Fst_2freq(data_ND_LASIDAD_EUR_EAS,col_NEA_EAS,col_NEA_EUR)$EAS_NEA_freq
data_ND_LASIDAD_EUR_EAS$Fst_NEA_EAS_LASIDAD<-Fst_2freq(data_ND_LASIDAD_EUR_EAS,col_NEA_EAS,col_NEA_LASIDAD)$EAS_NEA_freq
data_ND_LASIDAD_EUR_EAS$Fst_NEA_EUR_LASIDAD<-Fst_2freq(data_ND_LASIDAD_EUR_EAS,col_NEA_EUR,col_NEA_LASIDAD)$EUR_NEA_freq

##PBS
data_ND_LASIDAD_EUR_EAS$PBS_NEA <- (-log(1-data_ND_LASIDAD_EUR_EAS$Fst_NEA_EAS_LASIDAD)-log(1-data_ND_LASIDAD_EUR_EAS$Fst_NEA_EUR_LASIDAD)+log(1-data_ND_LASIDAD_EUR_EAS$Fst_NEA_EAS_EUR))/2

###Denisovan
##Fst
data_ND_LASIDAD_EUR_EAS$Fst_DEN_EAS_EUR<-Fst_2freq(data_ND_LASIDAD_EUR_EAS,col_DEN_EAS,col_DEN_EUR)$EAS_DEN_freq
data_ND_LASIDAD_EUR_EAS$Fst_DEN_EAS_LASIDAD<-Fst_2freq(data_ND_LASIDAD_EUR_EAS,col_DEN_EAS,col_DEN_LASIDAD)$EAS_DEN_freq
data_ND_LASIDAD_EUR_EAS$Fst_DEN_EUR_LASIDAD<-Fst_2freq(data_ND_LASIDAD_EUR_EAS,col_DEN_EUR,col_DEN_LASIDAD)$EUR_DEN_freq

##PBS
data_ND_LASIDAD_EUR_EAS$PBS_DEN <- (-log(1-data_ND_LASIDAD_EUR_EAS$Fst_DEN_EAS_LASIDAD)-log(1-data_ND_LASIDAD_EUR_EAS$Fst_DEN_EUR_LASIDAD)+log(1-data_ND_LASIDAD_EUR_EAS$Fst_DEN_EAS_EUR))/2


##Distribution PBS
data_w_N=filter(data_ND_LASIDAD_EUR_EAS,LASIDAD_NEA_freq>0 | EAS_NEA_freq>0 | EUR_NEA_freq>0)
data_w_D=filter(data_ND_LASIDAD_EUR_EAS,LASIDAD_DEN_freq>0 | EAS_DEN_freq>0 | EUR_DEN_freq>0)


##Neanderthal PBS distribution
quantile_95_NEA=quantile(data_w_N$PBS_NEA,0.95)
quantile_99_NEA=quantile(data_w_N$PBS_NEA,0.99)


#plot
ggplot(data_w_N,aes(x=PBS_NEA))+geom_histogram(aes(y = after_stat(count / sum(count))),binwidth=0.001)+
  xlab("PBS Neanderthal")+ylab("Frequency")+
  geom_vline(xintercept=quantile_99_NEA, color="red")+
  theme_light(base_size = 20)+
  coord_cartesian(xlim=c(-0.02,0.05))

##Denisova PBS distribution
quantile_95_DEN=quantile(data_w_D$PBS_DEN,0.95)
quantile_99_DEN=quantile(data_w_D$PBS_DEN,0.99)

##Plot distribution PBS
ggplot(data_w_D,aes(x=PBS_DEN))+geom_histogram(aes(y = after_stat(count / sum(count))),binwidth=0.0005)+
  xlab("PBS Denisovan")+ylab("Frequency")+ 
  geom_vline(xintercept=quantile_99_DEN, color="red")+
  theme_light(base_size = 20)+
  coord_cartesian(xlim=c(-0.005,0.05))



##group by chromosome/position
data_mutate_PBS_called = data_ND_LASIDAD_EUR_EAS %>%
  group_by(chrom, start) %>%
  summarize(
    mean_NEA_freq_LASIDAD=mean(LASIDAD_NEA_freq),
    mean_NEA_freq_EUR=mean(EUR_NEA_freq),
    mean_NEA_freq_EAS=mean(EAS_NEA_freq),
    mean_DEN_freq_LASIDAD=mean(LASIDAD_DEN_freq),
    mean_DEN_freq_EUR=mean(EUR_DEN_freq),
    mean_DEN_freq_EAS=mean(EAS_DEN_freq),
    mean_NEA_PBS = PBS_NEA, 
    mean_DEN_PBS = PBS_DEN,
    enriched_PBS_NEA = ifelse(called>0.5 &  LASIDAD_NEA_freq>EAS_NEA_freq & LASIDAD_NEA_freq>EUR_NEA_freq & LASIDAD_NEA_freq>lim_enriched_NEA & PBS_NEA>quantile_99_NEA,"1","0"), 
    enriched_PBS_DEN = ifelse(called>0.5 &  LASIDAD_DEN_freq>EAS_DEN_freq & LASIDAD_DEN_freq>EUR_DEN_freq & LASIDAD_DEN_freq>lim_enriched_DEN & PBS_DEN>quantile_99_DEN,"1","0") )%>%
  ungroup() %>%
  gather(archaic, freq, enriched_PBS_NEA:enriched_PBS_DEN) %>%
  mutate(Yaxis = ifelse(archaic == 'enriched_PBS_NEA', 1,2 )) 



data_mutate_PBS_called %>%
  mutate(chrom = factor(chrom, c(1:22) )) %>%
  
  ggplot() +
  
  # chromosomes positions
  geom_rect(data = chromosomes, aes(xmin = starts/1e6, xmax = ends/1e6, ymin = 0, ymax = length(Individuals) + 0.01), fill = 'white', color = 'white') +
  geom_text(data = chromosomes, aes(x = 0, y = (length(Individuals))/3, label = labels), hjust = 1) +
  geom_rect(data = gaps, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.02, ymax = length(Individuals) + 0.005), fill = 'white') +
  geom_rect(data = chromosomes, aes(xmin = starts/1e6, xmax = ends/1e6, ymin = 0, ymax = length(Individuals) + 0.005), fill = NA, color = 'white') +
  
  # Frequency in windows
  geom_rect(aes(xmin = start/1e6, xmax = (start+1000)/1e6, ymin = Yaxis-0.99, ymax = Yaxis-0.01, alpha=freq, fill=as.factor(Yaxis))) +
  
  ##add gaps
  geom_rect(data = gaps, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.02, ymax = length(Individuals) + 0.005), fill = 'white') +
  
  facet_grid(chrom~., switch = 'y') +
  
  scale_fill_manual(values = c("#56B4E9","#009E73"),breaks = c("1","2"), labels=c("Neanderthal","Denisovan"),name="Ancestry") +
  scale_alpha_manual(breaks=c("0","1"), labels=c("No","Yes"), values=c(0.1,1), name="Enriched (PBS>99%)")+
  theme_bw(base_size=15) +
  labs(x = 'genomic position (Mb)') +
  theme(strip.background = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_blank(), #element_text(angle = 180)
        legend.position = c(0.8, 0.2),
        panel.spacing.y=unit(0, "lines"))+
  labs(fill = "Ancestry")

#ggsave(paste0("LASIDAD_enriched_PBS_EAS_EUR_cutoff99_1kb_2sd_called05_EURphased_bluegreen.png"), width = 12, height = 5)
#ggsave("LASIDAD_enriched_PBS_EAS_EUR_cutoff99_1kb_2sd_called05_EURphased_bluegreen.pdf", width = 12, height = 5)
