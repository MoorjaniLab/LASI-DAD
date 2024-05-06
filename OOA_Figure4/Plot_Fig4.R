library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

options(dplyr.print_max = 1e9)
set.seed(1) 

pop_order = c('BantuSouthAfrica', 'Biaka', 'ESN', 'MSL', 'Yoruba', 'YRI', 'Mbuti', 'San', 'Mandenka', 'GWD', 'LWK', 'BantuKenya', 'ACB', 'ASW', 'CLM', 'PUR', 'MXL', 'PEL', 'Colombian', 'Karitiana', 'Maya', 'Pima', 'Surui', 'North', 'Central', 'West', 'South', 'East', 'North-East', 'Other', 'STU', 'BEB', 'GIH', 'ITU', 'PJL', 'Balochi', 'Brahui', 'Burusho', 'Hazara', 'Kalash', 'Makrani', 'Pathan', 'Sindhi', 'Uygur', 'JPT', 'KHV', 'CDX', 'CHB', 'CHS', 'Cambodian', 'Dai', 'Daur', 'Han', 'Hezhen', 'Japanese', 'Lahu', 'Miao', 'Mongolian', 'Naxi', 'NorthernHan', 'Oroqen', 'She', 'Tu', 'Tujia', 'Xibo', 'Yakut', 'Yi', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'Adygei', 'Basque', 'BergamoItalian', 'French', 'Orcadian', 'Russian', 'Sardinian', 'Tuscan', 'Bedouin', 'Druze', 'Mozabite', 'Palestinian', 'Bougainville', 'PapuanHighlands', 'PapuanSepik')


data = read.table('Data_Fig4.txt', header = T)

m1_all = lm(Emissions_human_phased~Archaic_sequence + missing, data )

INTERCEPT = m1_all$coefficients[[1]]
ARCHAIC = m1_all$coefficients[[2]]
MISSING = m1_all$coefficients[[3]]


CORRECTION = (0.45e-9 * 1000)

# ------------------------------------------------------------------------------------------------
# Correct for outgroup ancestry
corrected_data = data %>%
	mutate(value = (Emissions_human_phased - missing * MISSING - ARCHAIC * Archaic_sequence) / CORRECTION, facet = '2. Archaic ancestry') %>%
	select(pop, region, value, facet, dataset)


eruption_label = data.frame(x = 'Colombian', y = 75000, label = 'Estimated time of Toba Eruption (fossil record)')


corrected_data %>%
	mutate(color_scheme = ifelse(dataset == 'LASIDAD', 'India', as.character(region))) %>%
	mutate(pop = factor(pop, pop_order)) %>%

	ggplot() +
	geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 52781, ymax = 54918), fill = 'lightgrey') +
	geom_jitter(aes(x = pop, y = value, color = color_scheme), width = 0.25) +
	geom_hline(aes(yintercept = 74000), linetype = 'dashed') +
	geom_text(data = eruption_label, aes(x = x, y = y, label = label), hjust = 0) +
	theme_bw() +

	
	xlab('') +
	ylab('Minimum coalescence with outgroup (years)') +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
		  strip.text.x = element_blank(),
		  panel.grid.major.x = element_blank() ,
		  panel.spacing.x=unit(0.1, "lines"),
		  panel.spacing.y=unit(0.1, "lines"),
		  legend.position = 'bottom') +
	scale_color_manual('Regions',values = c(
				"AMERICA" = "#9A6324", 
				"SOUTH_ASIA" = "#469990",
				'India' = '#911eb4', 
				"EAST_ASIA" = "#000075",
				"EUROPE" = "#000000")) + 
	coord_cartesian(ylim = c(40000,80000))


ggsave('Figure4.png', width = 9, height = 7) 

