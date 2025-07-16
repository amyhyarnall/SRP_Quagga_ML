# Title: Step 4. SRP Plankton NMDS
# Author: Dr. Amy Yarnall
# Date Created: 2022-05-12
# Date last updated: 2025-07-15

# USEFUL LINKS
# https://stackoverflow.com/questions/25694515/convert-data-frame-to-a-matrix-with-column-1-of-df-as-rownames-of-matrix

# Load libraries
library(dplyr)
library(ggplot2)
library(vegan)
library(ggvegan) # this library needs a very new version of R
library(ggordiplots) # this library needs a very new version of R
library(gridExtra)
library(vegan3d)

# If you have trouble installing ggvegan or ggordiplots use:
# install.packages("remotes")
# remotes::install_github("gavinsimpson/ggvegan")
# remotes::install_github("jfq3/ggordiplots")

#### READ IN DATA AND WRANGLE ####
setwd('C:/Users/u4eewahy/Documents/BOR-SRP_Quagga_Habitat_Suitability/Ecosphere/GitHub data & code/Data')
phyto <- read.csv('ALL_SRP_Phytoplankton_data_20231030.csv')
zoo <- read.csv('ALL_SRP_Zooplankton_data_20231030.csv')

# Create column with waterbody ID
zoo$Waterbody <- ifelse(zoo$Station_ID %in% c("A-BHC","A-D","A-M","A-TMI"), "Apache Lake", 
                         ifelse(zoo$Station_ID %in% c("B-BC","B-D","B-M","B-MBL"), "Bartlett Reservoir", 
                                ifelse(zoo$Station_ID %in% c("C-BCF","C-BD","C-D","C-M"), "Canyon Lake",
                                       ifelse(zoo$Station_ID == "GR-D", "Granite Reef Diversion Reservoir",
                                              ifelse(zoo$Station_ID %in% c("R-CRS","R-D","R-M","R-WH"),"Theodore Roosevelt Lake","Saguaro Lake")))))
phyto$Waterbody <- ifelse(phyto$Station_ID %in% c("A-BHC","A-D","A-M","A-TMI"), "Apache Lake", 
                           ifelse(phyto$Station_ID %in% c("B-BC","B-D","B-M","B-MBL"), "Bartlett Reservoir", 
                                  ifelse(phyto$Station_ID %in% c("C-BCF","C-BD","C-D","C-M"), "Canyon Lake",
                                         ifelse(phyto$Station_ID == "GR-D", "Granite Reef Diversion Reservoir",
                                                ifelse(phyto$Station_ID %in% c("R-CRS","R-D","R-M","R-WH"),"Theodore Roosevelt Lake","Saguaro Lake")))))

# Create column with waterbody status
zoo$Status <- ifelse(zoo$Waterbody %in% c("Apache Lake","Canyon Lake","Saguaro Lake","Granite Reef Diversion Reservoir"), "Established","Negative")
phyto$Status <- ifelse(phyto$Waterbody %in% c("Apache Lake","Canyon Lake","Saguaro Lake","Granite Reef Diversion Reservoir"), "Established", "Negative")

wb_order <- c('Apache Lake','Bartlett Reservoir','Canyon Lake','Granite Reef Diversion Reservoir','Theodore Roosevelt Lake','Saguaro Lake')

# remove columns not needed from NMDS
drop_col <- c("Sample_Date","Sample_Time","Sample_Depth_m","Tally","NOTES", "Total_BV_um3.L", # to drop from zoo
              "Microcycstis","Comments","Tow_length_m","Net_radius_cm","Tow_vol_filtered_L",  # to drop from phyto
              "Total_sample_vol_mL","Aliquot_mL","Count_factor","Indiv_counted","Biomass_factor","Species_biomass_ugDW.L") 

zoo_reduced <- zoo[,!colnames(zoo) %in% drop_col]
phyto_reduced <- phyto[,!colnames(phyto) %in% drop_col]

# Create row ID for matrix
zoo_reduced$ID <- paste0(zoo_reduced$Station_ID, "_", zoo_reduced$Sample_Month, "-", zoo_reduced$Sample_Year)
phyto_reduced$ID <- paste0(phyto_reduced$Station_ID, "_", phyto_reduced$Sample_Month, "-", phyto_reduced$Sample_Year)

# Create new lowest taxa col for phyto and zoo
phyto_reduced$Lowest_taxon <- phyto_reduced$Genus
zoo_reduced$Lowest_taxon <- paste0(zoo_reduced$Genus, " ", zoo_reduced$Species)
zoo_reduced$Lowest_taxon <- gsub(" NA", "", zoo_reduced$Lowest_taxon)

# Combine different spellings of same species, assume unconfirmed identification ("cf.") are correct

# PHYTOPLANKTON
phyto_reduced$Lowest_taxon <- gsub("cf. ", "",phyto_reduced$Lowest_taxon) # remove unconfirmed designations
phyto_reduced$Lowest_taxon <- gsub(" spp.", " sp.",phyto_reduced$Lowest_taxon) # combine spp. and sp. designations
phyto_reduced$Lowest_taxon <- gsub("microscipicus", "microscopicus",phyto_reduced$Lowest_taxon) # fix misspelling 
phyto_reduced$Lowest_taxon <- gsub("Ulna", "ulna", phyto_reduced$Lowest_taxon) # fix misspelling 
phyto_reduced$Lowest_taxon <- gsub("vulgare", "vulgaris", phyto_reduced$Lowest_taxon) # alternate names for same species
phyto_reduced$Lowest_taxon <- gsub("minutulus", "minutus", phyto_reduced$Lowest_taxon) # alternate names for same species

# ZOOPLANKTON
zoo_reduced$Lowest_taxon <- gsub(" f. brevispinus", "",zoo_reduced$Lowest_taxon) # remove designations
zoo_reduced$Lowest_taxon <- gsub(" var. tecta", "",zoo_reduced$Lowest_taxon) # remove designations
zoo_reduced$Lowest_taxon <- gsub("trunctaum", "truncatum",zoo_reduced$Lowest_taxon) # fix misspelling 

# Put sample year into factor format
phyto_reduced$Sample_Year <- as.factor(phyto_reduced$Sample_Year)
zoo_reduced$Sample_Year <- as.factor(zoo_reduced$Sample_Year)

# Pull cols for NMDS
zoo_4sppmtrx <- zoo_reduced[,colnames(zoo_reduced) %in% c("ID", "Lowest_taxon", "Indiv_per_L")]
phyto_4sppmtrx  <- phyto_reduced[,colnames(phyto_reduced) %in% c("ID", "Lowest_taxon", "Density_cells.L")]

#### Create LOWEST TAXON plankton NMDS dataframes ####

# ZOOPLANKTON LOWEST TAXON
# create vectors for row and col names
zoospp_rownames <- unique(as.character(zoo_4sppmtrx$ID)) #Vector of rownames
zoospp_colnames <- unique(as.character(zoo_4sppmtrx$Lowest_taxon)) #Vector of colnames
# create matrix with cells = 0
zoo_spp_mtrx <- matrix(0, nrow = length(zoospp_rownames), ncol = length(zoospp_colnames))
row.names(zoo_spp_mtrx) <- zoospp_rownames
colnames(zoo_spp_mtrx) <- zoospp_colnames
# fill in matrix with the numbers by the index values from the df
zoo_spp_mtrx[as.matrix(zoo_4sppmtrx[c('ID','Lowest_taxon')])] <- as.double(zoo_4sppmtrx$Indiv_per_L)

# Add back in the columns for the meta data
zoo_spp_comm <- as.data.frame(cbind(ID = row.names(zoo_spp_mtrx), as.data.frame(zoo_spp_mtrx)))
zoo_spp_comm <- zoo_spp_comm[rowSums(zoo_spp_comm[,-1], na.rm = T)>0,] #remove missing sample
zoo_spp_comm <- zoo_spp_comm[,!names(zoo_spp_comm) == "veliger quagga"] # remove the quagga veligers from the dataset
zoo_spp_meta <- unique(zoo_reduced[, c("Waterbody","Status","Station_ID","Sample_Month","Sample_Year",'ID')])
zoo_spp_meta <- zoo_spp_meta[zoo_spp_meta$ID != 'C-BD_April-2023',] #remove missing sample
zoo_spp_NMDS <- merge(zoo_spp_meta, zoo_spp_comm, all = T, by = 'ID')

# PHYTOPLANKTON LOWEST TAXON
# create vectors for row and col names
phytospp_rownames <- unique(as.character(phyto_4sppmtrx$ID)) #Vector of rownames
phytospp_colnames <- unique(as.character(phyto_4sppmtrx$Lowest_taxon)) #Vector of colnames
# create matrix with cells = 0
phyto_spp_mtrx <- matrix(0, nrow = length(phytospp_rownames), ncol = length(phytospp_colnames))
row.names(phyto_spp_mtrx) <- phytospp_rownames
colnames(phyto_spp_mtrx) <- phytospp_colnames
# fill in matrix with the numbers by the index values from the df
phyto_spp_mtrx[as.matrix(phyto_4sppmtrx[c('ID','Lowest_taxon')])] <- as.double(phyto_4sppmtrx$Density_cells.L)

# Add back in the columns for the meta data
phyto_spp_comm <- as.data.frame(cbind(ID = row.names(phyto_spp_mtrx), as.data.frame(phyto_spp_mtrx)))
phyto_spp_meta <- unique(phyto_reduced[, c("Waterbody","Status","Station_ID","Sample_Month","Sample_Year",'ID')])
phyto_spp_NMDS <- merge(phyto_spp_meta, phyto_spp_comm, all = T, by = 'ID')

#### Transform  community data sets ####
# Step 1: remove any cols with sum of zero
# Step 2: fourth root data to minimize effects of zeros and extreme counts
# Step 3: Bray-Curtis distance within the NMDS

# zooplankon lowest taxon community - 4th root and Bray curtis
zoo_spp_comm <- zoo_spp_comm[,colSums(zoo_spp_comm[,-1], na.rm = T)>0]
zoo_spp_4rt <- (zoo_spp_comm[,-1])^(1/4)

# phytoplankon lowest taxon community - 4th root and Bray curtis
phyto_spp_comm <- phyto_spp_comm[,colSums(phyto_spp_comm[,-1], na.rm = T)>0]
phyto_spp_4rt <- (phyto_spp_comm[,-1])^(1/4)

#### NMDS: Phytoplankton lowest taxon level ####
set.seed(2022)
phyto_spp.nmds <- metaMDS(comm = phyto_spp_4rt, distance = 'bray', k = 2, try = 999,
                          engine = 'monoMDS', autotransform = F, noshare = T, stress = 1, 
                          wascores = T, expand = T, trace = F, plot = F)
phyto_spp.nmds$stress # <0.2
stressplot(phyto_spp.nmds)

phyto_spp.fit <- envfit(phyto_spp.nmds ~ Status, data = phyto_spp_NMDS, perm = 999); phyto_spp.fit
# Significant difference in centroids

# PLOT IT
phyto_spp.fort <- fortify(phyto_spp.nmds)
phyto_spp.fort.site <- subset(phyto_spp.fort, score == 'sites')
phyto_spp.fort.site$ID <- phyto_spp.fort.site$label
phyto_spp.fort.site <- merge(phyto_spp_meta, phyto_spp.fort.site, by = 'ID')

# basic 2D plot for reference only
phyto_spp_odriplot <- gg_ordiplot(phyto_spp.nmds, groups = phyto_spp_meta$Status, kind = "se", conf = 0.95, pt.size = 0, plot = T) 

phyto_ellipse <- phyto_spp_odriplot$df_ellipse
phytoEst_ellipse <- subset(phyto_ellipse, Group == "Established")
phytoNeg_ellipse <- subset(phyto_ellipse, Group == "Negative")

# Nicer 2D plot for final figure
phyto_spp_plot <- phyto_spp_odriplot$plot + theme_bw() +  xlim(-1.5,1.5) + ylim(-1.5,1.5) + 
  labs(x = 'NMDS1', y = 'NMDS2', col = "Status", shape = "Waterbody") +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', linewidth = 0.8) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', linewidth = 0.8) +
  geom_point(data = phyto_spp.fort.site, size = 2, 
             mapping = aes(x = NMDS1, y = NMDS2, col = Status, shape = Waterbody)) +
  geom_text(x = 1, y = 1.5,  # choose coordinate locations for stress metric text
            col = 'black', size = 4,
            aes(label = paste0("k = ",phyto_spp.nmds$ndim, ", Stress = ", round(phyto_spp.nmds$stress, digits = 3)))) +
  theme(legend.position = 'right') + 
  geom_polygon(data=phytoEst_ellipse, aes(x=x, y=y), fill ='lightcoral', col = "black", alpha = 0.4) +
  geom_polygon(data=phytoNeg_ellipse, aes(x=x, y=y), fill ='lightseagreen', col = "black", alpha = 0.4) +
  scale_color_manual(values = c('lightcoral','lightseagreen'))+
  scale_shape_manual(values = c(16,9,15,18,17,8))
phyto_spp_plot

#### NMDS: Zooplankton lowest taxon level ####
set.seed(2022)
zoo_spp.nmds <- metaMDS(comm = zoo_spp_4rt, distance = 'bray', k = 3, try = 999,
                        engine = 'monoMDS', autotransform = F, noshare = T, stress = 1, 
                        wascores = T, expand = T, trace = F, plot = F)
# NOTE: k = 3, stress = 0.25 when k = 2
zoo_spp.nmds$stress # <0.2 
stressplot(zoo_spp.nmds)
zoo_spp.fit <- envfit(zoo_spp.nmds ~ Status, data = zoo_spp_NMDS, perm = 999); zoo_spp.fit
# Significant difference in centroids

# PLOT IT
zoo_spp.fort <- fortify(zoo_spp.nmds)
zoo_spp.fort.site <- subset(zoo_spp.fort, score == 'sites')
zoo_spp.fort.site$ID <- zoo_spp.fort.site$label
zoo_spp.fort.site <- merge(zoo_spp_meta, zoo_spp.fort.site, by = 'ID')

# basic 2D plot for reference only
zoo_spp_odriplot <- gg_ordiplot(zoo_spp.nmds, groups = zoo_spp_meta$Status, kind = "se", conf = 0.95, pt.size = 0, plot = T) 

zoo_ellipse <- zoo_spp_odriplot$df_ellipse
zooEst_ellipse <- subset(zoo_ellipse, Group == "Established")
zooNeg_ellipse <- subset(zoo_ellipse, Group == "Negative")

# Nicer 2D plot for reference only
zoo_spp_plot <- zoo_spp_odriplot$plot + theme_bw() + 
  # CHOOSE the plot limits that are square and show all points
  xlim(-1.5,1.5) + ylim(-1.5,1.5) +
  labs(x = 'NMDS1', y = 'NMDS2', col = "Status", shape = "Waterbody") +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', linewidth = 0.8) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', linewidth = 0.8) +
  geom_point(data = zoo_spp.fort.site, size = 2,
             mapping = aes(x = NMDS1, y = NMDS2, col = Status, shape = Waterbody)) +
  geom_text(x = 1, y = 1.5, # choose coordinate locations for stress metric text
            col = 'black', size = 4,
            aes(label = paste0("k = ",zoo_spp.nmds$ndim, ", Stress = ",round(zoo_spp.nmds$stress, digits = 3)))) + 
  theme(legend.position = 'right') + 
  geom_polygon(data=zooEst_ellipse, aes(x=x, y=y), fill ='lightcoral', col = "black", alpha = 0.4) +
  geom_polygon(data=zooNeg_ellipse, aes(x=x, y=y), fill ='lightseagreen', col = "black", alpha = 0.4) +
  scale_color_manual(values = c('lightcoral','lightseagreen'))+
  scale_shape_manual(values = c(16,9,15,18,17,8)) 
zoo_spp_plot

# 3D plot - to be used in final figure
zoo_spp_odriplot3b <- with(zoo_spp_meta, ordiplot3d(zoo_spp.nmds,col=c('lightcoral','lightseagreen')[as.factor(Status)], pch=c(16,9,15,18,17,8)[as.factor(Waterbody)]))
with(zoo_spp_meta, ordiellipse(zoo_spp_odriplot3b, as.factor(Status), draw = "poly", col = c('lightcoral','lightseagreen'),
                               alpha = 0.4, label = F))

#### Figure 3. NMDS plot layout ####
# Function to separate legend for plotting
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}

layout1 <- rbind(c(1,1,2,2,3),
                 c(1,1,2,2,3))
grid.arrange(phyto_spp_plot + theme(legend.position = 'none') + ggtitle('A'),
             zoo_spp_plot +  theme(legend.position = 'none') + ggtitle('B'), # replace with 3d plot in final fig
             g_legend(zoo_spp_plot),
             layout_matrix = layout1)
# replace zoo 2d with 3d in layout
zoo_spp_odriplot3b <- with(zoo_spp_meta, ordiplot3d(zoo_spp.nmds,col=c('lightcoral','lightseagreen')[as.factor(Status)], pch=c(16,9,15,18,17,8)[as.factor(Waterbody)]))
with(zoo_spp_meta, ordiellipse(zoo_spp_odriplot3b, as.factor(Status), draw = "poly", col = c('lightcoral','lightseagreen'),
                               alpha = 0.4, label = F))

#### PERMANOVA and SIMPER analyses ####
# Run Permanova (adonis2 for Type II SS) for the lowest taxa level of phyto- and zoo-plankton
# Test effect of Status 
# Note: adonis wants the sample matrix and the metadata in different dataframes

# Use SIMPER to find out how dissimilar establish and negative communities are and which spp contribute most to dissimilarity
# Find whether spp of highest relative abundance contribute most to dissimilarity 

#****************************************
# Zooplankton lowest taxon level
# PERMANOVA - Status
adonis.zooltl1 <- adonis2(as.matrix(zoo_spp_4rt) ~ Status, data = zoo_spp_meta); adonis.zooltl1 

# SIMPER - Status
sim.zooltl.s <- with(zoo_spp_meta, simper(zoo_spp_4rt, group = Status, permutations = 100)) 

# What is the distribution of spp contributions to dissimilarities among statuses
par(mfrow = c(1,1))
hist(sim.zooltl.s$Established_Negative$average)

# It seems like most spp are contribution <3% so perhaps the few spp contributing >3% are more interesting for further analysis
# Which species are contributing (on average) > 3% to dissimilarity
# alphabetical order of top spp contributions
zoo_spp_diss <- sort(names(sim.zooltl.s$Established_Negative$average)[sim.zooltl.s$Established_Negative$average>0.03])
# What are all the unique species that contribute (on average) > 3% to dissimilarities btwn statuses
zoo_spp_diss
zoo_spp_diss_table <- data.frame(Mean_diss = sim.zooltl.s$Established_Negative$average, 
                                 SD = sim.zooltl.s$Established_Negative$sd)

# Beyond species contributing to dissimilarity, high (>1%) relative abundance of species may also be important
zoo_abun <- zoo_spp_NMDS %>% summarize_if(is.numeric, sum, na.rm=TRUE)
zoo_relabun <- (zoo_abun/rowSums(zoo_abun))*100

df_zoo_abun <- cbind(Species = rownames(t(zoo_relabun)), data.frame(t(zoo_relabun), row.names=NULL))
colnames(df_zoo_abun)[2]<-"Relative_abundance"
df_zoo_diss <- cbind(Species = rownames(zoo_spp_diss_table), data.frame(zoo_spp_diss_table*100, row.names=NULL))
colnames(df_zoo_diss)[2]<-"Mean_dissimilarity"
zoo_df <- merge(df_zoo_abun, df_zoo_diss)
zoo_df$Relative_abundance <- format(round(zoo_df$Relative_abundance, 4), nsmall = 4)
zoo_df$Mean_dissimilarity <- format(round(zoo_df$Mean_dissimilarity, 4), nsmall = 4)
zoo_df$SD <- format(round(zoo_df$SD, 4), nsmall = 4)

zoo_df <- zoo_df[with(zoo_df, order(Relative_abundance, Mean_dissimilarity, decreasing = T)), ]

par(mfrow = c(1,1))
hist(t(zoo_relabun), breaks = seq(0, 35, 1)) # the vast majority are <1%
zoo_relabun_1percent <- zoo_relabun[sapply(zoo_relabun, function(x) any(x >1))]
# these are the names of spp with >1% relative abundance
sort(names(zoo_relabun_1percent)) # note: includes spp contributing >3% dissimilarity are included 

zoo_relabun_table <- data.frame(Rel_abun = t(zoo_relabun_1percent))
zoo_table <- merge(zoo_relabun_table, zoo_spp_diss_table, by = "row.names", all = F)
colnames(zoo_table)[1] <- "Taxa"
zoo_table
#****************************************

#****************************************
# Phytoplankton lowest taxa level
# PERMANOVA - Status
adonis.phytoltl1 <- adonis2(as.matrix(phyto_spp_4rt) ~ Status, data = phyto_spp_meta); adonis.phytoltl1

# SIMPER - Status
sim.phytoltl.s <- with(phyto_spp_meta, simper(phyto_spp_4rt, group = Status, permutations = 100)) 

# What is the distribution of spp contributions to dissimilarities among statuses
par(mfrow = c(1,1))
hist(sim.phytoltl.s$Established_Negative$average)

data.frame(sim.phytoltl.s$Established_Negative$average)

# It seems like most spp are contribution <3% so perhaps the few spp contributing >3% are more interesting for further analysis
# Which species are contributing (on average) > 3% to dissimilarity
# alphabetical order of top spp contributions
phyto_spp_diss <- sort(names(sim.phytoltl.s$Established_Negative$average)[sim.phytoltl.s$Established_Negative$average>0.03])
# What are all the unique species that contribute (on average) > 3% to dissimilarities btwn statuses
phyto_spp_diss 
phyto_spp_diss_table <- data.frame(Mean_diss = sim.phytoltl.s$Established_Negative$average, 
                                   SD = sim.phytoltl.s$Established_Negative$sd)

# Beyond species contributing to dissimilarity, high (>1%) relative abundance of species may also be important
phyto_abun <- phyto_spp_NMDS %>% summarize_if(is.numeric, sum, na.rm=TRUE)
phyto_relabun <- (phyto_abun/rowSums(phyto_abun))*100

df_phyto_abun <- cbind(Species = rownames(t(phyto_relabun)), data.frame(t(phyto_relabun), row.names=NULL))
colnames(df_phyto_abun)[2]<-"Relative_abundance"
df_phyto_diss <- cbind(Species = rownames(phyto_spp_diss_table), data.frame(phyto_spp_diss_table*100, row.names=NULL))
colnames(df_phyto_diss)[2]<-"Mean_dissimilarity"
phyto_df <- merge(df_phyto_abun, df_phyto_diss)
phyto_df$Relative_abundance <- format(round(phyto_df$Relative_abundance, 4), nsmall = 4)
phyto_df$Mean_dissimilarity <- format(round(phyto_df$Mean_dissimilarity, 4), nsmall = 4)
phyto_df$SD <- format(round(phyto_df$SD, 4), nsmall = 4)

phyto_df <- phyto_df[with(phyto_df, order(Relative_abundance, Mean_dissimilarity, decreasing = T)), ]

par(mfrow = c(1,1))
hist(t(phyto_relabun), breaks = seq(0, 45, 1)) # the vast majority are <1%
phyto_relabun_1percent <- phyto_relabun[sapply(phyto_relabun, function(x) any(x >1))]
# these are the names of spp with >1% relative abundance
sort(names(phyto_relabun_1percent)) # note: includes spp contributing >3% dissimilarity are included 

phyto_relabun_table <- data.frame(Rel_abun = t(phyto_relabun_1percent))
phyto_table <- merge(phyto_relabun_table, phyto_spp_diss_table, by = "row.names", all = F)
colnames(phyto_table)[1] <- "Taxa"
phyto_table
#****************************************

#### Table S3, S4, S5 ####

# Table S3. Phytoplankton taxa observed listed in alphabetical order  
sort(unique(phyto_reduced$Lowest_taxon))

# Table S4. Zooplankton taxa observed listed in alphabetical order (note: excludes quagga veligers)
sort(unique(zoo_reduced$Lowest_taxon[zoo_reduced$Lowest_taxon != "veliger quagga"]))

# Table S5. Taxa that had >1% relative abundance with their respective phytoplankton and zooplankton communities. 
# Mean (SD) percent contribution to dissimilarity between established and negative 
# quagga population statuses from SIMPER analyses.

# Add divisions to table
phyto_div <- phyto_reduced %>%
  distinct(Lowest_taxon, Division, .keep_all = F)
colnames(phyto_div)[1]<- "Taxa"
zoo_div <- zoo_reduced %>%
  distinct(Lowest_taxon, Division, .keep_all = F)
colnames(zoo_div)[1]<- "Taxa"

phyto_table_div <- merge(phyto_table, phyto_div)
zoo_table_div <- merge(zoo_table, zoo_div)

# order by relative abundance
phyto_table_ord <- phyto_table_div[with(phyto_table_div, order(-Rel_abun)), ]
zoo_table_ord <- zoo_table_div[with(zoo_table_div, order(-Rel_abun)), ] 

# round numerics and order table cols
phyto_table_new <- data.frame(Taxa = phyto_table_ord$Taxa, Division = phyto_table_ord$Division, round(select(phyto_table_ord, -c(Taxa, Division)),2)) 
zoo_table_new <- data.frame(Taxa = zoo_table_ord$Taxa, round(select(zoo_table_ord, -c(Taxa, Division)),2)) 

# Table S5
phyto_table_new
zoo_table_new

# Write csvs with plankton names to use in step 4
# write.csv(sort(names(phyto_relabun_1percent)), 'ForGBM_phytoplankton_1percent_relativeabundance_20240124.csv', row.names = F)
# write.csv(sort(names(zoo_relabun_1percent)), 'ForGBM_zooplankton_1percent_relativeabundance_20240124.csv', row.names = F)
# ^ already available in data folder
