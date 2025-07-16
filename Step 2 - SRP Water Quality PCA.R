# Title: Step 2. Water quality PCA
# Author: Dr. Amy Yarnall
# Date Created: 2022-05-12
# Date last updated: 2025-07-15

# USEFUL LINKS
# https://www.datacamp.com/community/tutorials/pca-analysis-r
# https://tem11010.github.io/Plotting-PCAs/

# Load libraries
library(FactoMineR)
library(ggplot2)
library(cowplot)
library(dplyr)

#### READ IN DATA AND WRANGLE ####
# Set working directory and load datasets
setwd('C:/Users/u4eewahy/Documents/BOR-SRP_Quagga_Habitat_Suitability/Ecosphere/GitHub data & code/Data')
chl <- read.csv('ALL_SRP_Chlorophyll_data_20231030.csv')
mp <- read.csv('ALL_SRP_Multiprobe_data_20231030.csv')
wq <- read.csv('ALL_SRP_WQ_report_data_20231030.csv')

# Remove metadata and excluded var cols from datasets
chl_reduced <- chl[,!colnames(chl) %in% c("Lab_Number", "Vol_filtered_mL", "Sample_Date", "Chlb_mg.m3","Chlc_mg.m3")] 
mp_reduced <- mp[,!colnames(mp) %in% c("Sample_Date","Sample_Time","Site","mmHg","DO_percent","SPC_uS.cm","pH","pH_mV","DEP_m")]
wq_reduced <- wq[,!colnames(wq) %in% c("Lab_Number","Sampled_by","Sample_Date","Sample_Time","Date_received",    
                                       "pH","Sum_TDS_ppm","TSS_ppm","Na_ppm","Mg_ppm","CO3_ppm",          
                                       "HCO3_ppm","Cl_ppm","Total_Dissolved_P","NH3_N","NO3_NO2_asN",
                                       "As","Be","Cd","Cr","Cu","Mn","Mo","Ni","Pb","Zn","Al","Fe","Se")]

## Consistent sample IDs and formatting needed to merge dataset

## Chlorophyll - only surface ('-S') samples
# create separate column for location in water column
chl_reduced$Location <- ifelse(grepl('-S', chl_reduced$Station_ID), "Surface", "Bottom") 
colnames(chl_reduced)[1] <- 'Station_ID_location'
# remove '-S' tag from sample ID
chl_reduced$Station_ID <- gsub('.{2}$', "", chl_reduced$Station_ID_location)  

## Water quality - both bottom ('-B') and surface ('-S') samples
wq_reduced$Location <- ifelse(grepl('-S', wq_reduced$Station_ID), "Surface", "Bottom") 
colnames(wq_reduced)[1] <- 'Station_ID_location'
# remove '-S' and '-B' tags from sample ID
wq_reduced$Station_ID <- gsub('.{2}$', "", wq_reduced$Station_ID_location) 

## Multiparameter probe - samples throughout water column
# Take  water column mean
mp_reduced$Deg_C <- as.numeric(mp_reduced$Deg_C) # covert character to numeric
mp_reduced_watercolmeans <- mp_reduced %>% group_by(Sample_Year, Sample_Month, Station_ID) %>%
  summarise(colmu_Deg_C = mean(as.numeric(Deg_C, na.rm = T)),
            colmu_DO_mg.L = mean(DO_mg.L, na.rm = T))
# Create separate column for location in water column
mp_reduced_watercolmeans$Location <- rep("Full water column", length(mp_reduced_watercolmeans$Station_ID))


## Take  means of SI converted data at station level

chl_reduced <- chl_reduced %>% mutate(across(ends_with("_mg.m3"), ~ .x / 1000, .names = "{.col%>%sub('_mg.m3', '_mg.L', .)}"))
chl_reduced$Pheophytin_mg.L <- chl_reduced$Pheophytin/1000
chl_station <- chl_reduced %>% group_by(Station_ID) %>%
  summarise(mu_Chla_mg.L = mean(Chla_mg.L, na.rm = T),
            mu_Pheophytin_mg.L = mean(Pheophytin_mg.L, na.rm = T), 
            mu_Chla_Pheo_664.665a = mean(Chla_Pheo_664.665a, na.rm = T))

wq_reduced <- wq_reduced %>% mutate(across(ends_with("_ppm"), ~ .x * 0.9988590004, .names = "{.col%>%sub('_ppm', '_mg.L', .)}"))
wq_reduced <- wq_reduced %>% mutate(across(all_of(c("Total_P","Total_N")), ~ .x * 0.9988590004, .names = "{paste0(.col, '_mg.L')}"))
wq_station <- wq_reduced %>% group_by(Station_ID) %>%
  summarise(mu_TDS_mg.L = mean(TDS_mg.L, na.rm = T), 
            mu_K_mg.L = mean(K_mg.L, na.rm = T), 
            mu_Ca_mg.L = mean(Ca_mg.L, na.rm = T),
            mu_Alk_mg.L = mean(Alk_mg.L, na.rm = T), 
            mu_SO4_mg.L = mean(SO4_mg.L, na.rm = T), 
            mu_SiO2_mg.L = mean(SiO2_mg.L, na.rm = T), 
            mu_F_mg.L = mean(F_mg.L, na.rm = T), 
            mu_TP_mg.L = mean(Total_P, na.rm = T), 
            mu_TN_mg.L = mean(Total_N, na.rm = T))

mp_station <- mp_reduced_watercolmeans %>% group_by(Station_ID) %>%
  summarise(mu_Deg_C = mean(colmu_Deg_C, na.rm = T), 
            mu_DO_mg.L = mean(colmu_DO_mg.L, na.rm = T))

## MERGE DATASETS
srp.wq <- merge(chl_station,merge(wq_station, mp_station, all = T), all = T)


# Create column with waterbody ID
srp.wq$Waterbody <- ifelse(srp.wq$Station_ID %in% c("A-BHC","A-D","A-M","A-TMI"), "Apache Lake", 
                            ifelse(srp.wq$Station_ID %in% c("B-BC","B-D","B-M","B-MBL"), "Bartlett Reservoir", 
                                   ifelse(srp.wq$Station_ID %in% c("C-BCF","C-BD","C-D","C-M"), "Canyon Lake",
                                          ifelse(srp.wq$Station_ID == "GR-D", "Granite Reef Diversion Reservoir",
                                                 ifelse(srp.wq$Station_ID %in% c("R-CRS","R-D","R-M","R-WH"),"Theodore Roosevelt Lake","Saguaro Lake")))))
# Create column with waterbody status
srp.wq$Status <- ifelse(srp.wq$Waterbody %in% c("Apache Lake","Canyon Lake","Saguaro Lake","Granite Reef Diversion Reservoir"), "Established", "Negative")
# Create dataframe only with meta-data
srp.meta <- select(srp.wq, Waterbody, Status, Station_ID)
# Create dataframe only with data
srp.data <- srp.wq[,!colnames(srp.wq) %in% colnames(srp.meta)]
# Recombine into new DF with new column order
srp.means.pca <- cbind(srp.meta, srp.data)
# Create unique row names
row.names(srp.means.pca) <- srp.means.pca$Station_ID

#### FactoMineR PCA Station Means (all data averaged per Station ID) ####
# How much data are we inputting into the PCA?
dim(srp.means.pca[,!colnames(srp.means.pca) %in% colnames(srp.meta)])
# 20 stations, 14 WQ variables 

# Run the PCA
srp.PCA <- PCA(srp.means.pca[,!colnames(srp.means.pca) %in% colnames(srp.meta)], #run the PCA on the data (exclude metadata cols)
               scale.unit = T, # scale each variable to a z-score 
               graph = F) # Don't generate the graph, this will be done later with ggplot 

summary(srp.PCA) # examine the results 
srp.PCA$eig[1:2,] # eigenvalues and % variance accounted for in 2D

# Extract pc scores for first two components and add to meta-dataframe
srp.meta$PC1 <- srp.PCA$ind$coord[, 1] # indexing the first column
srp.meta$PC2 <- srp.PCA$ind$coord[, 2]  # indexing the second column

# Extract the data for the variable contributions to each of the pc axes
pca.vars <- srp.PCA$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
rownames(pca.vars) <- 1:length(pca.vars$vars)

# Subset to the variables that load on a dimension >0.4
# By convention this means they are contributing enough to the Dimensions 1 or 2 to merit further exploration 
pca.vars.sub <- pca.vars[,colnames(pca.vars) %in% c('vars', 'Dim.1', 'Dim.2')]
pca.vars.sub <- pca.vars.sub[(abs(pca.vars.sub$Dim.1) > 0.4 | abs(pca.vars.sub$Dim.2) > 0.4), ]
dim(pca.vars.sub)[1] # Number of variables that load on a dim.1 or dim.2  >0.4

# Create more readable labels for the plot from the var names
unique(pca.vars.sub$vars)

pca.vars.sub$labels <- gsub("mu_", '', pca.vars.sub$vars) # remove mean designation
pca.vars.sub$labels <- gsub("_mg.L", '', pca.vars.sub$labels) # remove unit
pca.vars.sub$labels <- gsub("Deg_C", 'Temperature', pca.vars.sub$labels) # replace this label
pca.vars.sub$labels <- gsub("Alk", 'Alkalinity', pca.vars.sub$labels) # replace this label
pca.vars.sub$labels <- gsub("Chla_Pheo_664.665a", 'Chl a/Pheophytin', pca.vars.sub$labels) # replace this label
pca.vars.sub$labels <- gsub("_", ' ', pca.vars.sub$labels) # replace remaining _ with spaces
pca.vars.sub$labels <- gsub("Chl", 'Chl ', pca.vars.sub$labels) # add space after Chl

unique(pca.vars.sub$labels)

#### Figure 2. PCA Plot ####
# This function will be used in the creation of the PCA plots
# By convention, the variable contribution plot has a circle around the variables that has a radius of 1. 
# Here's a function to make one.
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
circ <- circleFun(c(0,0),2,npoints = 500)
options(ggrepel.max.overlaps = Inf)

p1 <- ggplot() + 
  geom_point(data = srp.meta, aes(x = PC1, y = PC2, shape = Waterbody, color = Status), size = 2) +  
  stat_ellipse(data = srp.meta, geom = "polygon", aes(x = PC1, y = PC2, color = Status), 
               alpha = 0, lwd = 1.2, show.legend = FALSE, level = 0.95) +
  geom_hline(yintercept = 0, lty = 2) + 
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = paste0('PC1 (',round(srp.PCA$eig[1,2],2),'%)'), # percentage of variance accounted for in Dim 1
       y = paste0('PC2 (',round(srp.PCA$eig[2,2],2),'%)'), # percentage of variance accounted for in Dim 2
       color = 'Status', shape = 'Waterbody') + 
  scale_color_manual(values = c('lightcoral','lightseagreen')) +
  scale_shape_manual(values = c(16,9,15,18,17,8)) +
  theme_minimal() + 
  theme(panel.grid = element_blank(), panel.border = element_rect(fill= "transparent"), 
        legend.position = "right")

p.vars <- ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(data = pca.vars.sub, aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 1, alpha = 0.5) + 
  ggrepel::geom_text_repel(data = pca.vars.sub, 
                           aes(x = Dim.1, y =  Dim.2, 
                               label = labels), size = 4, point.padding = 0.5) +
  labs(x = paste0('PC1 (',round(srp.PCA$eig[1,2],2),'%)'), 
       y = paste0('PC2 (',round(srp.PCA$eig[2,2],2),'%)')) + 
  coord_equal() +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

prow1 <- plot_grid(p1+theme(legend.position = "none"), 
                   p.vars, rel_heights = c(1, 2), rel_widths = c(1.5, 1.5), labels = c("A", "B"))

# Function to separate legend for plotting
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}
# Plot with ellipses
plot_grid(prow1, g_legend(p1), nrow = 1, rel_widths = c(2, 0.5), rel_heights = c(1, 1))
# Manually adjust dimensions of plot area to get proper layout

#### Table S2. PCA loadings table ####
# Table of loading values sorted by absolute value of PC1 then PC2
loadings <- pca.vars.sub[,c(4,1,2)]
names(loadings) <- c("Water quality parameter", "PC1", "PC2")
loadings_ord <- loadings[with(loadings, order(-abs(PC1), -abs(PC2))), ]
loadings_ord$PC1 <- round(loadings_ord$PC1, 2)
loadings_ord$PC2 <- round(loadings_ord$PC2, 2)
loadings_ord

# All 14 WQ vars make it on to step 4
# Write csvs with wq var names to use in step 4
# write.csv(data.frame(vars = pca.vars$vars), 'ForGBM_WQvars_0.4loadingPCA_20250528.csv', row.names = F)
# ^ already available in data folder