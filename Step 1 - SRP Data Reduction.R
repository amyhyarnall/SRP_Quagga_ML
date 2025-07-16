# Title: Step 1. SRP data reduction
# Author: Dr. Amy Yarnall
# Date Created: 2022-05-12
# Date last updated: 2025-07-15

# Load libraries
library(dplyr)
library(ggplot2)
library(caret)

#### APPENDIX S1: Supplemental analyses ####
# READ IN DATA
# Set working directory and load datasets
setwd('C:/Users/u4eewahy/Documents/BOR-SRP_Quagga_Habitat_Suitability/Ecosphere/GitHub data & code/Data')

# Read in the latest versions of the 
# Water quality parameters - for PCA
chl <- read.csv('ALL_SRP_Chlorophyll_data_20231030.csv', fileEncoding = 'UTF-8-BOM')
mp <- read.csv('ALL_SRP_Multiprobe_data_20231030.csv', fileEncoding = 'UTF-8-BOM')
wq <- read.csv('ALL_SRP_WQ_report_data_20231030.csv', fileEncoding = 'UTF-8-BOM')

# Quagga mussel veliger count data
quagga <- read.csv('ALL_SRP_Mussel_data_20250527.csv', fileEncoding = 'UTF-8-BOM')

# Reduce each of these data sets only to the needed columns
# Remove metadata cols from datasets
chl_reduced <- chl[,!colnames(chl) %in% c('Lab_Number', 'Vol_filtered_mL')] 
mp_reduced <- mp[,!colnames(mp) == 'Sample_Time']
wq_reduced <- wq[,!colnames(wq) %in% c('Lab_Number', 'Sampled_by', 'Sample_Time', 'Date_received')]
quagga_reduced <- quagga[,!colnames(quagga) %in% c('Lab_number', 'Location')]

# Remove veliger samples not yet analyzed
quagga_reduced$Veligers_L <- as.numeric(quagga_reduced$Veligers_L) # converts column to numeric, makes "Not yet analyzed" in to NAs
quagga_reduced <- quagga_reduced[!is.na(quagga_reduced$Veligers_L),] # removes NAs from dataset

#### FIGURE S1: VELIGER COUNT DATA ####
mussel_station_sum <- quagga_reduced %>% group_by(Water_body, Station_ID, Sample_Year) %>%
  summarise(yr_sum = sum(Veligers_L, na.rm = T), yr_mean = mean(Veligers_L, na.rm = T), yr_sd = sd(Veligers_L, na.rm = T))
mussel_waterbody_sum <- mussel_station_sum %>% group_by(Water_body, Sample_Year) %>%
  summarise(mean_st_velcount = mean(yr_sum), sd_st_velcount = sd(yr_sum), mean_veliger = mean(yr_mean), sd_veliger = sd(yr_mean))

# VISUALIZE MUSSEL VELIGER DATA
# reorder reservoirs
quagga_reduced$Water_body <- factor(quagga_reduced$Water_body, 
                                    levels = c("Apache Lake", "Canyon Lake", "Saguaro Lake", "Granite Reef Diversion Dam", 
                                               "Bartlett Reservoir", "Theodore Roosevelt Lake"))
levels(quagga_reduced$Water_body)[4] <- "Granite Reef Diversion Reservoir" #change waterbody name
# Plot
ggplot(quagga_reduced) +
  geom_boxplot(aes(x = Water_body, y = log10(Veligers_L), col = Water_body)) +
  labs(x = NULL, y = expression(Log[10]~veligers/L)) + theme_classic() +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "#00BFC4")) +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1), legend.position = "none") + facet_wrap(~Sample_Year)

#### Table S1: DATA REDUCTION MEASURES 1 and 2 ####

# 1. All measures of pH were excluded from analysis due to:
# a. probe calibration errors
# b. time elapsed from sample collection to measurement in lab

# 2. Some related groups of parameters were reduced to single parameters, because they were redundant 
# (e.g., multiple units of measure) or were partial components of another parameter 
# a. DO (mg/L) & DO (%) - redundant
# b. TDS (ppm) & Sum TDS (ppm) - redundant
# c. Total dissolved P = major partial component of TP
# d. NH3 as N (ppm) & NO3/NO2 as N (ppm) = major partial components of TN

#### Table S1: DATA REDUCTION MEASURE 3 ####

# 3. Exclude any remaining ex situ variables with > or = 20% of samples at their lower method detection limit (MDL)

## Chlorophyll data
# What % of samples are at the lower MDL (0.35)
length(chl_reduced$Chla_mg.m3[chl_reduced$Chla_mg.m3 <= 0.35])/length(chl_reduced$Chla_mg.m3) # 0%
length(chl_reduced$Chlb_mg.m3[chl_reduced$Chlb_mg.m3 <= 0.35])/length(chl_reduced$Chlb_mg.m3) # 68% - exclude
length(chl_reduced$Chlc_mg.m3[chl_reduced$Chlc_mg.m3 <= 0.35])/length(chl_reduced$Chlc_mg.m3) # 22% - exclude
length(chl_reduced$Pheophytin[chl_reduced$Pheophytin <= 0.35])/length(chl_reduced$Pheophytin) # 10%

# Per Chemistry lab - Remove April 2023 no3/no2 samples > total N; these are likely sampling errors
# Excluded all NH3 and no3/no2 samples > total N
wq_reduced$NO3_NO2_asN[wq_reduced$NO3_NO2_asN > wq_reduced$Total_N] <- NA
wq_reduced$NH3_N[wq_reduced$NH3_N > wq_reduced$Total_N] <- NA

# What % of WQ samples are at their lower method detection limits (L-MDL)
length(wq_reduced$TDS_ppm[wq_reduced$TDS_ppm == 2.00])/length(wq_reduced$TDS_ppm) # 0.7%
length(wq_reduced$K_ppm[wq_reduced$K_ppm %in% c(1.56, 0.29)])/length(wq_reduced$K_ppm) # 0%
length(wq_reduced$Ca_ppm[wq_reduced$Ca_ppm %in% c(1.26, 2.03)])/length(wq_reduced$Ca_ppm) # 0%
length(wq_reduced$Alk_ppm[wq_reduced$Alk_ppm == 1.4])/length(wq_reduced$Alk_ppm) # 0.7%
length(wq_reduced$SO4_ppm[wq_reduced$SO4_ppm %in% c(0.2, 0.27, 0.15)])/length(wq_reduced$SO4_ppm) # 0%
length(wq_reduced$SiO2_ppm[wq_reduced$SiO2_ppm %in% c(0.59, 1.09)])/length(wq_reduced$SiO2_ppm) # 0%
length(wq_reduced$F_ppm[wq_reduced$F_ppm %in% c(0.04, 0.01)])/length(wq_reduced$F_ppm) # 0%
length(wq_reduced$Total_P[wq_reduced$Total_P %in% c(0.005,0.003)])/length(wq_reduced$Total_P) # 1%
length(wq_reduced$Total_N[wq_reduced$Total_N %in% c(0.04,0.006)])/length(wq_reduced$Total_N) # 0%
length(wq_reduced$TSS_ppm[wq_reduced$TSS_ppm == 2.00])/length(wq_reduced$TSS_ppm) # 32% - exclude
length(wq_reduced$Cl_ppm[wq_reduced$Cl_ppm %in% c(0.06, 0.31)])/length(wq_reduced$Cl_ppm) # 0%
length(wq_reduced$CO3_ppm[wq_reduced$CO3_ppm == 1.4])/length(wq_reduced$CO3_ppm) # 84% - exclude
length(wq_reduced$Total_Dissolved_P[wq_reduced$Total_Dissolved_P %in% c(0.005,0.003)])/length(wq_reduced$Total_Dissolved_P) # 18%
length(wq_reduced$NH3_N[wq_reduced$NH3_N %in% c(0.00600, 0.00200)])/length(wq_reduced$NH3_N) # 9%
length(wq_reduced$NO3_NO2_asN[wq_reduced$NO3_NO2_asN %in% c(0.00800, 0.00400)])/length(wq_reduced$NO3_NO2_asN) # 21% - exclude

# Metals - exclude all
length(wq_reduced$Al[wq_reduced$Al %in% c(2.86, 181)])/length(wq_reduced$Al) # 87%
length(wq_reduced$As[wq_reduced$As %in% c(30.7, 8.12)])/length(wq_reduced$As) # 97%
length(wq_reduced$Be[wq_reduced$Be  %in% c(1.2, 3.36)])/length(wq_reduced$Be) # 100%
length(wq_reduced$Cd[wq_reduced$Cd  %in% c(1.6, 7.04)])/length(wq_reduced$Cd) # 100%
length(wq_reduced$Cr[wq_reduced$Cr  %in% c(1.67, 5.44)])/length(wq_reduced$Cr) # 100%
length(wq_reduced$Cu[wq_reduced$Cu  %in% c(2.02, 0.57)])/length(wq_reduced$Cu) # 41%
length(wq_reduced$Fe[wq_reduced$Fe %in% c(2.22, 1.63)])/length(wq_reduced$Fe) # 20%
length(wq_reduced$Mn[wq_reduced$Mn %in% c(14.2, 2.57)])/length(wq_reduced$Mn) # 82%
length(wq_reduced$Mo[wq_reduced$Mo %in% c(2.8, 4.06)])/length(wq_reduced$Mo) # 98%
length(wq_reduced$Ni[wq_reduced$Ni %in% c(8.16, 7.66)])/length(wq_reduced$Ni) # 100%
length(wq_reduced$Pb[wq_reduced$Pb %in% c(8.0, 29.0)])/length(wq_reduced$Pb) # 100%
length(wq_reduced$Se[wq_reduced$Se %in% c(27.8, 23.5)])/length(wq_reduced$Se) # 91%
length(wq_reduced$Zn[wq_reduced$Zn %in% c(19.2, 544)])/length(wq_reduced$Zn) # 96%

#### Table S1: DATA REDUCTION MEASURE 4 ####

# 4. Find any highly correlated variables 

# First each dataset needs consistent sample IDs and formatting to merge and compare
# Also convert cols to SI units

## Chlorophyll - only surface ('-S') samples
# create separate column for location in water column
chl_reduced$Location <- ifelse(grepl('-S', chl_reduced$Station_ID), "Surface", "Bottom") 
colnames(chl_reduced)[1] <- 'Station_ID_location'
# remove '-S' tag from sample ID
chl_reduced$Station_ID <- gsub('.{2}$', "", chl_reduced$Station_ID_location)  
# When are these samples taken?
chl_reduced$Sample_mo_yr <- paste0(chl_reduced$Sample_Month, "_", chl_reduced$Sample_Year)
# Create columns that convert original units to SI units
chl_reduced <- chl_reduced %>% mutate(across(ends_with("_mg.m3"), ~ .x / 1000, .names = "{.col%>%sub('_mg.m3', '_mg.L', .)}"))
chl_reduced$Pheophytin_mg.L <- chl_reduced$Pheophytin/1000
  # ratio col does not need converting because units cancelled 
# Remove cols with original units
chl_reduced <- chl_reduced[!colnames(chl_reduced) %in% c('Chla_mg.m3','Chlb_mg.m3','Chlc_mg.m3','Pheophytin')]

## Water quality - both bottom ('-B') and surface ('-S') samples
wq_reduced$Location <- ifelse(grepl('-S', wq_reduced$Station_ID), "Surface", "Bottom") 
colnames(wq_reduced)[1] <- 'Station_ID_location'
# remove '-S' and '-B' tags from sample ID
wq_reduced$Station_ID <- gsub('.{2}$', "", wq_reduced$Station_ID_location) 
# Sample Month_Year
wq_reduced$Sample_mo_yr <- paste0(wq_reduced$Sample_Month, "_", wq_reduced$Sample_Year)
# Create columns that convert original units to SI units - only for remaining considered vars
# Starting unit = ppm or ppb
wq_reduced <- wq_reduced %>% mutate(across(ends_with("_ppm"), ~ .x * 0.9988590004, .names = "{.col%>%sub('_ppm', '_mg.L', .)}"))
wq_reduced <- wq_reduced %>% mutate(across(all_of(c("Total_P","Total_N","Total_Dissolved_P","NH3_N","NO3_NO2_asN")), ~ .x * 0.9988590004, .names = "{paste0(.col, '_mg.L')}"))
wq_reduced <- wq_reduced %>% mutate(across(all_of(c('Al','As','Be','Cd','Cr','Cu','Fe','Mn','Mo','Ni','Pb','Se','Zn')), ~ .x * 0.009988590004, .names = "{paste0(.col, '_mg.L')}"))
## Remove cols with original units
wq_reduced <- wq_reduced[!colnames(wq_reduced) %in% c("Total_P","Total_N","Total_Dissolved_P","NH3_N","NO3_NO2_asN",'Al','As','Be','Cd','Cr','Cu','Fe','Mn','Mo','Ni','Pb','Se','Zn')]
wq_reduced <- wq_reduced %>% select(-matches("_ppm"))

## Multiprobe - samples throughout water column
# First pool depth measurements to get a water column mean - excluding vars from previous data reduction steps
mp_reduced_watercolmeans <- mp_reduced %>% group_by(Sample_Year, Sample_Month, Station_ID) %>%
  summarise(n_DEP = length(DEP_m), max_DEP_m = max(DEP_m, na.rm = T), # n measurements and max depth
            colmu_Deg_C = mean(as.numeric(Deg_C, na.rm = T)),
            min_mmHg = min(mmHg, na.rm = T), max_mmHg = max(mmHg, na.rm = T), # min and max pressure readings
            colmu_DO_percent = mean(DO_percent, na.rm = T),
            colmu_DO_mg.L = mean(DO_mg.L, na.rm = T), 
            colmu_SPC_uS.cm = mean(SPC_uS.cm, na.rm = T),
            colmu_pH = mean(pH, na.rm = T), 
            colmu_pH_mV = mean(pH_mV, na.rm = T))

# Create separate column for location in water column
mp_reduced_watercolmeans$Location <- rep("Full water column", length(mp_reduced_watercolmeans$Station_ID))
# Sample Month_Year
mp_reduced_watercolmeans$Sample_mo_yr <- paste0(mp_reduced_watercolmeans$Sample_Month, "_", mp_reduced_watercolmeans$Sample_Year)


# Take mean of each variable at level of STATION ID (PER YEAR-MONTH) and merge to single dataframe
chl_yrmo <- chl_reduced %>% group_by(Station_ID,Sample_Year,Sample_Month,Sample_mo_yr) %>%
  summarise(mu_Chla_mg.L = mean(Chla_mg.L, na.rm = T),
            mu_Chlb_mg.L = mean(Chlb_mg.L, na.rm = T),
            mu_Chlc_mg.L = mean(Chlc_mg.L, na.rm = T),
            mu_Pheophytin_mg.L = mean(Pheophytin_mg.L, na.rm = T), 
            mu_Chla_Pheo_664.665a = mean(Chla_Pheo_664.665a, na.rm = T))

wq_yrmo <- wq_reduced %>% group_by(Station_ID,Sample_Year,Sample_Month,Sample_mo_yr) %>%
  summarise(mu_EC_uS.cm = mean(EC_uS.cm, na.rm = T), 
            mu_wq_pH = mean(pH, na.rm = T), 
            mu_TDS_mg.L = mean(TDS_mg.L, na.rm = T), 
            mu_Sum_TDS_mg.L = mean(Sum_TDS_mg.L, na.rm = T), 
            mu_TSS_mg.L = mean(TSS_mg.L, na.rm = T),
            mu_Na_mg.L = mean(Na_mg.L, na.rm = T), 
            mu_K_mg.L = mean(K_mg.L, na.rm = T), 
            mu_Ca_mg.L = mean(Ca_mg.L, na.rm = T),
            mu_Mg_mg.L = mean(Mg_mg.L, na.rm = T), 
            mu_CO3_mg.L = mean(CO3_mg.L, na.rm = T), 
            mu_HCO3_mg.L = mean(HCO3_mg.L, na.rm = T), 
            mu_Alk_mg.L = mean(Alk_mg.L, na.rm = T), 
            mu_Cl_mg.L = mean(Cl_mg.L, na.rm = T), 
            mu_SO4_mg.L = mean(SO4_mg.L, na.rm = T), 
            mu_SiO2_mg.L = mean(SiO2_mg.L, na.rm = T), 
            mu_F_mg.L = mean(F_mg.L, na.rm = T), 
            mu_TP_mg.L = mean(Total_P_mg.L, na.rm = T), 
            mu_Total_Dissolved_P_mg.L = mean(Total_Dissolved_P_mg.L, na.rm = T), 
            mu_NH3_N_mg.L = mean(NH3_N_mg.L, na.rm = T), 
            mu_NO3_NO2_asN_mg.L = mean(NO3_NO2_asN_mg.L, na.rm = T), 
            mu_TN_mg.L = mean(Total_N_mg.L, na.rm = T), 
            mu_As_mg.L = mean(As_mg.L, na.rm = T), 
            mu_Be_mg.L = mean(Be_mg.L, na.rm = T), 
            mu_Cd_mg.L = mean(Cd_mg.L, na.rm = T), 
            mu_Cr_mg.L = mean(Cr_mg.L, na.rm = T), 
            mu_Cu_mg.L = mean(Cu_mg.L, na.rm = T), 
            mu_Mn_mg.L = mean(Mn_mg.L, na.rm = T), 
            mu_Mo_mg.L = mean(Mo_mg.L, na.rm = T), 
            mu_Ni_mg.L = mean(Ni_mg.L, na.rm = T), 
            mu_Pb_mg.L = mean(Pb_mg.L, na.rm = T), 
            mu_Zn_mg.L = mean(Zn_mg.L, na.rm = T), 
            mu_Al_mg.L = mean(Al_mg.L, na.rm = T), 
            mu_Fe_mg.L = mean(Fe_mg.L, na.rm = T), 
            mu_Se_mg.L = mean(Se_mg.L, na.rm = T))

mp_yrmo <- mp_reduced_watercolmeans %>% group_by(Station_ID,Sample_Year,Sample_Month,Sample_mo_yr) %>%
  summarise(mu_n_DEP = mean(n_DEP, na.rm = T), mu_max_DEP_m = mean(max_DEP_m, na.rm = T),
            mu_Deg_C = mean(colmu_Deg_C, na.rm = T), 
            mu_mmHg = mean(max_mmHg, na.rm = T), 
            mu_DO_percent = mean(colmu_DO_percent, na.rm = T), 
            mu_DO_mg.L = mean(colmu_DO_mg.L, na.rm = T), 
            mu_SPC_uS.cm = mean(colmu_SPC_uS.cm, na.rm = T), 
            mu_mp_pH = mean(colmu_pH, na.rm = T), 
            mu_mp_pH_mV = mean(colmu_pH_mV, na.rm = T))
# Check that all allowable sample IDs remain

table(chl_yrmo$Station_ID, chl_yrmo$Sample_mo_yr) 
# 2021: Apr, Jul; 2022: Jan, Apr, Jul, Oct; 2023: Jan, Apr, Jul
# known missing sample: A-D on April_2021

table(wq_yrmo$Station_ID, wq_yrmo$Sample_mo_yr) 
# 2021: Apr, Jul; 2022: Jan, Apr, Jul, Oct; 2023: Jan, Apr, Jul
# known missing sample: A-D on April_2021

table(mp_yrmo$Station_ID, mp_yrmo$Sample_mo_yr) 
# Sampled more often but contains all nine sample periods from other datasets
# known missing station samples throughout
# 2021: Feb, Mar, Apr, May, Jul, Aug, Sep, Dec
# 2022: Jan, Mar, Apr, Jul, Aug, Sep, Oct, Nov, Dec
# 2023: Jan, Feb, Mar, Apr, May, Jun, Jul, Aug

## MERGE DATASETS
Station_ID.yrmo <- merge(chl_yrmo,
                       merge(wq_yrmo, mp_yrmo, all = T), all = T)
# only keep the months where all chl, wq, and mp samples were taken
shared_months <- Reduce(intersect, list(unique(chl_yrmo$Sample_mo_yr),
                                        unique(wq_yrmo$Sample_mo_yr),
                                        unique(mp_yrmo$Sample_mo_yr)))
Station_ID.yrmo_reduced <- Station_ID.yrmo[Station_ID.yrmo$Sample_mo_yr %in% shared_months,]


## Explore correlations among remaining variables
# It may be helpful to reduce the dimensionality of this dataset PRIOR to PCA
# Including a lot of highly correlated variables in PCA can inflate the similarity between points

vars <- select(Station_ID.yrmo_reduced, -c(Sample_Year, Sample_Month, Sample_mo_yr, Station_ID))
corr <- cor(na.omit(vars))
highlyCorr <- findCorrelation(corr, cutoff = .8)
corrvars <- names(vars[,highlyCorr]); corrvars # names of vars that have high correlations to others

# Which variables are correlated to each other - find all pairs
corr_table <- data.frame(corr)
corr_table[corr_table == 1] <- NA
corr_table[abs(corr_table) > 0.8] <- 'high'
corr_high <- which(corr_table=='high', arr.in= TRUE)
highcorr_table <- cbind.data.frame(var1 = colnames(corr_table)[ corr_high[, 1] ],
                                   var2 = rownames(corr_table)[ corr_high[, 2] ])
# Shows individual variable and all correlates
highcorr_groups <- highcorr_table %>% group_by(variable = var1) %>% summarise(correlates = paste(var2, collapse = ", "))

# Which variables are correlated to the most others?
highcorr_counts <- data.frame(table(highcorr_table$var1))
hist(highcorr_counts$Freq) 
colnames(highcorr_counts)[1] <- 'variable'
highcorr_count <- merge(highcorr_counts, highcorr_groups, all = F); highcorr_count

# Correlated groups of vars should be reduced to one representative

#### FIGURE S2:  CORRELATED VARIABLE GROUPS ####
corr_fig <- Station_ID.yrmo_reduced

colClean <- function(x) {colnames(x) <- gsub("mu_", "", colnames(x)) # remove mu tag
colnames(x) <- gsub("_uS.cm", "", colnames(x)) # remove unit
colnames(x) <- gsub("_mg.L", "", colnames(x)) # remove unit
colnames(x) <- gsub(" ","_",  colnames(x))} #change spaces to underscores - more r friendly

colnames(corr_fig) <- colClean(corr_fig)
corr_fig$Status <- as.factor(ifelse(corr_fig$Station_ID %in% c("B-BC","B-D","B-M","B-MBL","R-CRS","R-D","R-M","R-WH"),"Negative", "Established"))

featurePlot(x = corr_fig[,names(corr_fig) %in% c('TDS','SPC','Na','Cl')], 
            y = corr_fig$Status, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)), 
            layout = c(5,1))

featurePlot(x = corr_fig[,names(corr_fig) %in% c("Alk","Mg","HCO3")], 
            y = corr_fig$Status, 
            plot = "box", 
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)), 
            layout = c(5,1))

featurePlot(x = corr_fig[,names(corr_fig) %in% c("SO4", "Mg")], 
            y = corr_fig$Status, 
            plot = "box", 
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)), 
            layout = c(5,1))

## For step 2 and 4 KEEP TDS, ALKALINITY, AND SO4, Discard other correlated vars

# Write csv for use in step 4
# write.csv(Station_ID.yrmo_reduced, "SRP_StationID_YearMonth_SI_20250528.csv", row.names = F)
# ^ already available in data folder