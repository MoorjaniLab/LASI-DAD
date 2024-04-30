#####
## Plots Figure 2 - HBD // IBD
#####
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggnewscale)

###
## Metadata
###
#Region
col_regions<-c("North"="#f58231","North-East"="#4363d8","East"="#42d4f4","Central"="#e6194B","West"="#ffe119","South"="#3cb44b","Other"="#a9a9a9",
                  "SAS"="#469990","EUR"="#800000","EAS"="#000000","AFR"="#808000","LASI-DAD"="#911eb4", "LASI-DAD bootstrap"="#911eb4")

col_regions_LASIDAD<-c("North"="#f58231","North-East"="#4363d8","East"="#42d4f4","Central"="#e6194B","West"="#ffe119","South"="#3cb44b","Other"="#a9a9a9",
                       "SAS"="grey","EUR"="grey","EAS"="grey","AFR"="grey","LASI-DAD"="grey", "LASI-DAD bootstrap"="grey")

order_regions<-c("AFR","EAS","EUR","SAS","North","North-East","East","Central","West","South","Other")

label_regions<-c("North"="LASI-DAD North","North-East"="LASI-DAD North-East","East"="LASI-DAD East","Central"="LASI-DAD Central","West"="LASI-DAD West","South"="LASI-DAD South","Other"="LASI-DAD Other",
                 "SAS"="1000G SAS","EUR"="1000G EUR","EAS"="1000G EAS","AFR"="1000G AFR")
col_regions_cousins<-c("North"="#f58231","North-East"="#4363d8","East"="#42d4f4","Central"="#e6194B","West"="#ffe119","South"="#3cb44b","Other"="#a9a9a9",
                       "SAS"="#469990","EUR"="#800000","EAS"="#000000","AFR"="#808000","LASI-DAD"="#911eb4", "LASI-DAD bootstrap"="#911eb4",
                       "1st degree"="#ffd8b1","2nd degree"="#dcbeff","3rd degree"='#fabed4',"4th degree"='#aaffc3' ,"5th degree"='#fffac8')

#State 
order_states<-c("Jammu_&_Kashmir","Punjab","Haryana","Delhi","Rajasthan","Himachal_Pradesh",
                "Uttarakhand","Uttar_Pradesh","Madhya_Pradesh",
                "Gujarat","Maharashtra",
                "Telangana","Karnataka","Tamil_Nadu","Kerala","Andhra_Pradesh",
                "Bihar","West_Bengal","Odisha","Jharkhand",
                "Assam","Meghalaya")
label_state<-c("Jammu_&_Kashmir"="Jammu & Kashmir","Punjab"="Punjab","Uttarakhand"="Uttranchal","Haryana"="Haryana","Delhi"="Delhi",
               "Rajasthan"="Rajasthan","Uttar_Pradesh"="Uttar Pradesh",
               "Bihar"="Bihar","West_Bengal"="West Bengal","Assam"="Assam","Odisha"="Orissa","Gujarat"="Gujarat","Madhya_Pradesh"="Madhya Pradesh",
               "Maharashtra"="Maharashtra","Telangana"="Telangana","Karnataka"="Karnataka","Tamil_Nadu"="Tamil Nadu","Kerala"="Kerala",
               "Meghalaya"="Meghalaya","Himachal_Pradesh"="Himachal Pradesh","Jharkhand"="Jharkhand","Andhra_Pradesh"="Andhra Pradesh")
col_state<-c("Jammu_&_Kashmir"="#8C564B","Punjab"="#FFBB78","Uttarakhand"="#F7B6D2","Haryana"="#9467BD","Delhi"="#2CA02C",
             "Rajasthan"="#98DF8A","Uttar_Pradesh"="#C49C94","Bihar"="#FF7F0E","West_Bengal"="#C7C7C7","Assam"="#1F77B4","Odisha"="#AEC7E8",
             "Gujarat"="#D62728","Madhya_Pradesh"="#BCBD22","Maharashtra"="#17BECF","Telangana"="#C5B0D5","Karnataka"="#E377C2","Tamil_Nadu"="#FF9896","Kerala"="#7F7F7F",
             "Meghalaya"="blue","Himachal_Pradesh"="red","Jharkhand"="green","Andhra_Pradesh"="gold")
#Bootstrap
dashtype=c("bootstrap"="dotted","full sample"="solid")

#Kth degree cousin postions
pos_x_cousin<-c("1st degree"=847.752,"2nd degree"=221.94,"3rd degree"=52.98,"4th degree"=13.25,"5th degree"=3.31)
text_cousins=c("1st degree","2nd degree","3rd degree","4th degree","5th degree" )


###
## HBD
###

##Data
file_hbd_LASIDAD="LASI_DAD2620_HBD_sup2cM_sup8cM_sup20cM.txt"
hbd_LASIDAD=read_table(file_hbd_LASIDAD,col_names = TRUE)


file_hbd_1000G="1000G_nochild_HBD_sup2cM_sup8cM_sup20cM.txt"
hbd_1000G=read_table(file_hbd_1000G,col_names = TRUE)

merged_dataset <- rbind(hbd_LASIDAD,hbd_1000G)

## Proportions
#over 8cM
ggplot()+
  geom_boxplot(data=filter(hbd_LASIDAD,State!='abroad') ,aes(State,(nb_sup2-nb_sup20)/(nb_sup2),fill=State, col=Region))+
  scale_fill_manual(labels=label_state,values=col_state,limits=order_states)+
  scale_color_manual(values=col_regions,limits=c("North","Central","West","South","East","North-East"))+
  scale_x_discrete(labels=label_state,limits=order_states)+
  ylab("Proportion of HBD segments <8cM")+xlab("")+
  theme_light(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#over 8cM
ggplot()+
  geom_boxplot(data=filter(hbd_LASIDAD,State!='abroad') ,aes(State,(nb_sup2-nb_sup20)/(nb_sup2),fill=State, col=Region))+
  scale_fill_manual(labels=label_state,values=col_state,limits=order_states)+
  scale_color_manual(values=col_regions,limits=c("North","Central","West","South","East","North-East"))+
  scale_x_discrete(labels=label_state,limits=order_states)+
  ylab("Proportion of HBD segments <20cM")+xlab("")+
  theme_light(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

###Bar plot, amount per individual
#>8cM
ggplot(data = filter(merged_dataset,Region!='AMR' & Region != "Other" & Region!="North-East"))+
  geom_bar(aes(x=reorder(`#ID`,-sum_sup2),y=sum_sup2,fill=Region,alpha="below_8"),stat = "identity")+
  geom_bar(aes(x=reorder(`#ID`,-sum_sup2),y=sum_sup8,alpha="above_8"),fill="black",stat = "identity")+
  scale_alpha_manual(values=c("below_8"=1, "above_8"=0.8), labels=c("above_8"="> 8cM"),breaks=c("above_8"),
                     name="Lengths of homozygous segments")+
  ylab("Total homozygous sequence (cM)")+xlab("")+
  facet_grid(~fct_relevel(Region,'South','North','Central','East','West','SAS','EAS','EUR','AFR'),
             scales = "free", space='free',
             switch="both")+
  scale_fill_manual(labels=label_regions,values=col_regions_LASIDAD,limits=order_regions)+
  scale_y_log10()+
  theme_minimal(base_size = 30)+
  guides(fill=FALSE)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.grid.major = element_blank(),
        legend.position="bottom",
        legend.key = element_rect(fill = "white", colour = "black"),
        strip.background =element_blank(),
        strip.text = element_text(colour = 'black'))

#>20cM
ggplot(data = filter(merged_dataset,Region!='AMR' & Region != "Other" & Region!="North-East"))+
  geom_bar(aes(x=reorder(`#ID`,-sum_sup2),y=sum_sup2,fill=Region,alpha="below_8"),stat = "identity")+
  geom_bar(aes(x=reorder(`#ID`,-sum_sup2),y=sum_sup20,alpha="above_8"),fill="black",stat = "identity")+
  scale_alpha_manual(values=c("below_8"=1, "above_8"=0.8), labels=c("above_8"="> 20cM"),breaks=c("above_8"),
                     name="Lengths of homozygous segments")+
  ylab("Total homozygous sequence (cM)")+xlab("")+
  facet_grid(~fct_relevel(Region,'South','North','Central','East','West','SAS','EAS','EUR','AFR'),
             scales = "free", space='free',
             switch="both")+
  scale_fill_manual(labels=label_regions,values=col_regions_LASIDAD,limits=order_regions)+
  scale_y_log10()+
  theme_minimal(base_size = 30)+
  guides(fill=FALSE)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.grid.major = element_blank(),
        legend.position="bottom",
        legend.key = element_rect(fill = "white", colour = "black"),
        strip.background =element_blank(),
        strip.text = element_text(colour = 'black'))

###
## IBD
###

##Data
file_ibd_close_LASIDAD="LASI_DAD2620_IBD_closest_sup2.txt"
ibd_close_LASIDAD=read_table(file_ibd_close_LASIDAD,col_names = T)

file_ibd_close_1000G="1000G_IBD_closest_sup2.txt"
ibd_close_1000G=read_table(file_ibd_close_1000G,col_names = T)

merged_dataset_ibd_close <- rbind(ibd_close_LASIDAD,filter(ibd_close_1000G,Region1=="AFR" | Region1=="EUR" | Region1 == "EAS" | Region1=="SAS"))

##Boxplot IBD shared with closest relative
ggplot()+
  geom_boxplot(data = merged_dataset_ibd_close, aes(x=Region1,y=sum_sup2,col=Region1),notch=TRUE)+
  scale_color_manual(labels=label_regions,values=col_regions_cousins,limits=order_regions)+
  scale_x_discrete(limits=order_regions,labels=label_regions)+
  new_scale_color()+
  geom_hline(aes(yintercept=pos_x_cousin, col=text_cousins),linetype = "dashed",linewidth=1)+
  scale_color_manual(values=col_regions_cousins,guide="none")+
  geom_text(aes(label=text_cousins,y=exp(log(pos_x_cousin)+0.11),color=text_cousins),x=11,size=4)+
  scale_y_log10()+
  xlab("")+ylab("Sum of IBD segments shared with closest related individual")+
  theme_light(base_size = 20)+
  theme(legend.title=element_blank())


## Bootstrap
file_ibd_close_n500_LASIDAD="closest_for_500ind_500times_percolumn.txt"
ibd_close_n500_LASIDAD=read_table(file_ibd_close_n500_LASIDAD,col_names = F)

# Find info 95% distribution for 50%
row_50p=ibd_close_n500_LASIDAD[250,2:1001]
lim025=quantile(unlist(row_50p),0.025,names=F)
med=quantile(unlist(row_50p),0.5,names=F)
lim975=quantile(unlist(row_50p),0.975,names=F)

#95%-median for every row--> create new dataset
#init
percentile_025=c()
percentile_50=c()
percentile_79=c()
percentile_975=c()

for (n in 1:500){
  row=ibd_close_n500_LASIDAD[n,2:1001]
  percentile_025[n] <- quantile(unlist(row),0.025,names=F)
  percentile_50[n] <- quantile(unlist(row),0.5,names=F)
  percentile_79[n] <- quantile(unlist(row),0.79,names=F)
  percentile_975[n] <- quantile(unlist(row),0.975,names=F)
}


##Plot cumulative 1000G + LASI-DAD and Boostrap
ggplot() +
  geom_line(data=filter(ibd_close_1000G,Region1!="AMR"), aes(x = sum_sup2,col=Region1,y = (1 - ..y.. )*100,linetype="full sample"), stat='ecdf',linewidth=1)+
  geom_line(data=ibd_close_LASIDAD, aes(x = sum_sup2,y = (1 - ..y.. )*100,col="LASI-DAD",linetype="full sample"), stat='ecdf',linewidth=1)+
  geom_line(aes(x=percentile_025,y=ibd_close_n500_LASIDAD$X1*100,col="LASI-DAD",linetype="bootstrap"),linewidth=0.5)+
  geom_line(aes(x=percentile_50,y=ibd_close_n500_LASIDAD$X1*100,col="LASI-DAD",linetype="bootstrap"),linewidth=1)+
  geom_line(aes(x=percentile_975,y=ibd_close_n500_LASIDAD$X1*100,col="LASI-DAD",linetype="bootstrap"),linewidth=0.5)+
  scale_color_manual(labels=label_regions,values=col_regions)+
  scale_linetype_manual(values=dashtype)+
  scale_x_log10()+
  ylab("Individuals (%)")+xlab("IBD sharing cutoff (CM) to closest related individual")+
  new_scale_color()+
  geom_vline(aes(xintercept=pos_x_cousin, col=text_cousins),linetype = "dashed",linewidth=1.5)+
  scale_color_manual(values=col_regions_cousins,guide = 'none')+
  geom_text(aes(label=text_cousins,x=pos_x_cousin+1,color=text_cousins),y=104,size=5)+
  theme_light(base_size = 20)+
  theme(legend.title=element_blank())

##Different SSU
file_ibd_closediffssuid_LASIDAD="LASI_DAD2620_IBD_closest_diffssu_sup2.txt"
ibd_closediffssuid_LASIDAD=read_table(file_ibd_closediffssuid_LASIDAD,col_names = T)

ggplot() +
  scale_color_manual(labels=label_regions,values=col_regions_cousins,limits=order_regions)+
  geom_line(data=ibd_closediffssuid_LASIDAD, aes(x = sum_sup2,col=Region1,y = (1 - ..y.. )*100,linetype = "in different ssu"), stat='ecdf',linewidth=1)+
  geom_line(data=ibd_close_LASIDAD, aes(x = sum_sup2,col=Region1,y = (1 - ..y.. )*100), stat='ecdf',linewidth=1)+
  geom_line(data=filter(ibd_close_1000G,Region1!="AMR"), aes(x = sum_sup2,col=Region1,y = (1 - ..y.. )*100), stat='ecdf',linewidth=1)+
  scale_x_log10()+
  scale_linetype_manual(breaks=c("in different ssu"), values = c("dashed"))+
  ylab("Individual (%)")+xlab("IBD sharing cutoff (cM) to closest related individual")+
  new_scale_color()+
  geom_vline(aes(xintercept=pos_x_cousin, col=text_cousins),linetype = "dashed",linewidth=1)+
  scale_color_manual(values=col_regions_cousins,guide = 'none')+
  geom_text(aes(label=text_cousins,x=pos_x_cousin+1,color=text_cousins),y=104,size=5)+  theme_light(base_size = 20)+
  theme(legend.title=element_blank())