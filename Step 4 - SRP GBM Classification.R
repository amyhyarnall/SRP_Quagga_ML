# Title: Step 5. SRP Classification GBM (Gradient boosted machine)
# Author: Dr. Amy Yarnall
# Date Created: 2022-07-13
# Date Last Updated: 2025-07-15

# Load libraries
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation
library(ggplot2) # for plotting
# library(ggpattern) # for adding pattern to geom_bar
library(gridExtra) # for plotting
library(grid) # for plotting
# library(rpart) # for tree
# library(rattle) # for plotting the tree
library(caret) # for gradient boosted modeling
library(gbm) # for examining variable importance
# library(MLmetrics) 
library(pROC) # for ROC curves
# library(stringr) # for file name string manipulation
library(pdp) # for partial dependence plots
# library(psych)

# Helpful links - If you want a primer on gradient boosted decision trees, watch these videos 
# https://www.youtube.com/watch?v=3CC4N4z3GJc
# https://www.youtube.com/watch?v=2xudPOBz-vs
# https://www.youtube.com/watch?v=jxuNLH5dXCs
# https://www.youtube.com/watch?v=StWY5QWMXCw

# Helpful links - tutorials for caret:
# https://topepo.github.io/caret/index.html # This one is most thorough
# https://quantdev.ssri.psu.edu/sites/qdev/files/09_EnsembleMethods_2017_1127.html
# https://www.projectpro.io/recipes/apply-gradient-boosting-for-classification-r
# https://www.youtube.com/watch?v=GZafoSe5DvI
# https://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
# https://www.machinelearningplus.com/machine-learning/caret-package/

####### Section 1: Load data ####
setwd('C:/Users/u4eewahy/Documents/BOR-SRP_Quagga_Habitat_Suitability/Ecosphere/GitHub data & code/Data')

# Read in dataframes that contain the names of the variables that were important from the PCA and NMDSs
wq_vars <- read.csv('ForGBM_WQvars_0.4loadingPCA_20250528.csv', fileEncoding = 'UTF-8-BOM')
phyto_vars <- read.csv('ForGBM_phytoplankton_1percent_relativeabundance_20240124.csv', fileEncoding = 'UTF-8-BOM')
zoo_vars <- read.csv('ForGBM_zooplankton_1percent_relativeabundance_20240124.csv', fileEncoding = 'UTF-8-BOM')

# Read dataframes that contain the data for the variables names to be put into the classification tree
srp.wqyrmo <- read.csv('SRP_StationID_YearMonth_SI_20250528.csv', fileEncoding = 'UTF-8-BOM')
srp.phyto <- read.csv('ALL_SRP_Phytoplankton_data_20231030.csv', fileEncoding = 'UTF-8-BOM')
srp.zoo <- read.csv('ALL_SRP_Zooplankton_data_20231030.csv', fileEncoding = 'UTF-8-BOM')

# read in response variable (mussel) data
srp.mussel <- read.csv('ALL_SRP_Mussel_data_20250527.csv', fileEncoding = 'UTF-8-BOM')
colnames(srp.mussel)[3] <- "Waterbody" # Change colname from Water_body 

####### Section 2: Summarize WQ data and plankton data to the year-month level ####
# Create metadata df
meta_yrmo <- unique(srp.mussel[names(srp.mussel) %in% c('Waterbody', 'Station_ID', 'Sample_Year','Sample_Month')])
meta_yrmo$Sample_mo_yr <- paste0(meta_yrmo$Sample_Month,"_",meta_yrmo$Sample_Year)

# Established and Negative statuses only
meta_yrmo$Status <- ifelse(meta_yrmo$Waterbody %in% c("Apache Lake","Canyon Lake","Saguaro Lake","Granite Reef Diversion Dam"), "Established", "Negative")

# update Granite reef name
meta_yrmo$Waterbody <- ifelse(meta_yrmo$Waterbody %in% c("Granite Reef Diversion Dam"), "Granite Reef Diversion Reservoir", meta_yrmo$Waterbody)
srp.mussel$Waterbody <- ifelse(srp.mussel$Waterbody %in% c("Granite Reef Diversion Dam"), "Granite Reef Diversion Reservoir", srp.mussel$Waterbody)

# Summarize mussel columns needed for GBM
mussel_yrmo <- srp.mussel %>% group_by(Waterbody,Station_ID,Sample_Year,Sample_Month) %>%
  summarise(sum_veligers_L = sum(as.numeric(Veligers_L), na.rm = T), 
            mu_veligers_L = mean(as.numeric(Veligers_L), na.rm = T))
mussel_yrmo$Sample_mo_yr <- paste0(mussel_yrmo$Sample_Month,"_",mussel_yrmo$Sample_Year)
mussel_yrmo$ID <- paste0(mussel_yrmo$Station_ID,"_",mussel_yrmo$Sample_mo_yr)

# Pull out Water Quality columns from data files based upon names listed in the wq_vars dataframes
wq_reduced_yrmo <- srp.wqyrmo[names(srp.wqyrmo) %in% c('Station_ID','Sample_Year','Sample_Month','Sample_mo_yr',wq_vars$vars)]
wq_reduced_yrmo$ID <- paste0(wq_reduced_yrmo$Station_ID,"_",wq_reduced_yrmo$Sample_mo_yr)

# Pull out Plankton counts from data files based upon names listed in the _vars dataframes

# phytoplankton
# Combine different spellings of same species, assume unconfirmed identification ("cf.") are correct
srp.phyto$Genus <- gsub("cf. ", "",srp.phyto$Genus) # remove unconfirmed designations
srp.phyto$Genus <- gsub("microscipicus", "microscopicus",srp.phyto$Genus) # fix misspelling 

phyto <- srp.phyto[srp.phyto$Genus %in% phyto_vars$x, ]
phyto$Lowest_taxon <- phyto$Genus
phyto$Sample_mo_yr <- paste0(phyto$Sample_Month,"_",phyto$Sample_Year)

# zooplankton
srp.zoo$Lowest_taxon <- paste0(srp.zoo$Genus, " ", srp.zoo$Species)
srp.zoo$Lowest_taxon <- gsub(" NA", "", srp.zoo$Lowest_taxon)
# Combine different spellings of same species
srp.zoo$Lowest_taxon <- gsub(" var. tecta", "", srp.zoo$Lowest_taxon) # remove designations

zoo <- srp.zoo[srp.zoo$Lowest_taxon %in% zoo_vars$x, ]
zoo$Sample_mo_yr <- paste0(zoo$Sample_Month,"_",zoo$Sample_Year)

# phytoplankton - year month summary
phyto_yrmosum <- phyto %>% group_by(Station_ID, Sample_Year, Sample_Month, Sample_mo_yr, Lowest_taxon) %>%
  summarise(mu_Density_cells.L = mean(Density_cells.L))

# zooplankton - year month summary
zoo_yrmosum <- zoo %>% group_by(Station_ID, Sample_Year, Sample_Month, Sample_mo_yr, Lowest_taxon) %>%
  summarise(mu_Indiv_per_L = mean(Indiv_per_L))

# Reformat each plankton species to be columns

# PHYTOPLANKTON - Year Month
phyto_yrmosum$ID <- paste0(phyto_yrmosum$Station_ID,"_",phyto_yrmosum$Sample_mo_yr)
phyto_4mtrx2 <- phyto_yrmosum[,colnames(phyto_yrmosum) %in% c("ID","Lowest_taxon","mu_Density_cells.L")]
# create vectors for row and col names
phyto_rownames2 <- unique(as.character(phyto_yrmosum$ID)) #Vector of rownames
phyto_colnames2 <- unique(as.character(phyto_yrmosum$Lowest_taxon)) #Vector of colnames
# create matrix with cells = 0
phyto_mtrx2 <- matrix(0, nrow = length(phyto_rownames2), ncol = length(phyto_colnames2))
row.names(phyto_mtrx2) <- phyto_rownames2
colnames(phyto_mtrx2) <- phyto_colnames2
# fill in matrix with the numbers by the index values from the df
phyto_mtrx2[as.matrix(phyto_4mtrx2[c('ID','Lowest_taxon')])] <- as.double(phyto_4mtrx2$mu_Density_cells.L)
# make a df from the matrix
phyto_comm2 <- as.data.frame(cbind(ID = row.names(phyto_mtrx2), as.data.frame(phyto_mtrx2)))
# add units to colnames
colnames(phyto_comm2)[-1] <- paste(colnames(phyto_comm2)[-1],"no.L",sep="_")

# ZOOPLANKTON - Year Month
zoo_yrmosum$ID <- paste0(zoo_yrmosum$Station_ID,"_",zoo_yrmosum$Sample_mo_yr)
zoo_4mtrx2 <- zoo_yrmosum[,colnames(zoo_yrmosum) %in% c("ID","Lowest_taxon","mu_Indiv_per_L")]
# create vectors for row and col names
zoo_rownames2 <- unique(as.character(zoo_yrmosum$ID)) #Vector of rownames
zoo_colnames2 <- unique(as.character(zoo_yrmosum$Lowest_taxon)) #Vector of colnames
# create matrix with cells = 0
zoo_mtrx2 <- matrix(0, nrow = length(zoo_rownames2), ncol = length(zoo_colnames2))
row.names(zoo_mtrx2) <- zoo_rownames2
colnames(zoo_mtrx2) <- zoo_colnames2
# fill in matrix with the numbers by the index values from the df
zoo_mtrx2[as.matrix(zoo_4mtrx2[c('ID','Lowest_taxon')])] <- as.double(zoo_4mtrx2$mu_Indiv_per_L)
# make a df from the matrix
zoo_comm2 <- as.data.frame(cbind(ID = row.names(zoo_mtrx2), as.data.frame(zoo_mtrx2)))
# add units to colnames
colnames(zoo_comm2)[-1] <- paste(colnames(zoo_comm2)[-1],"no.L",sep="_")

# ONLY WQ DATA AND SAMPLE MONTHS
gbm_wqinput_yrmo <- merge(meta_yrmo, 
                           merge(mussel_yrmo, wq_reduced_yrmo, all = T), all = T)
gbm_wqinput_yrmo <- gbm_wqinput_yrmo[gbm_wqinput_yrmo$Sample_mo_yr %in% unique(wq_reduced_yrmo$Sample_mo_yr),]

# ONLY PLANKTON DATA AND shared SAMPLE MONTHS
plankton_sharedmonths <- Reduce(intersect, list(unique(phyto_yrmosum$Sample_mo_yr),
                                                unique(zoo_yrmosum$Sample_mo_yr)))
gbm_planktoninput_yrmo <- merge(meta_yrmo,
                                 merge(mussel_yrmo, 
                                       merge(phyto_comm2, zoo_comm2, all = T), all = T), all = T)
gbm_planktoninput_yrmo <- gbm_planktoninput_yrmo[gbm_planktoninput_yrmo$Sample_mo_yr %in% plankton_sharedmonths,]

####### Section 3: Final formatting of each dataset for Station Year-Month gbms ####
# Make the IDs the rownames
rownames(gbm_wqinput_yrmo) <- gbm_wqinput_yrmo$ID 
rownames(gbm_planktoninput_yrmo) <- gbm_planktoninput_yrmo$ID 

# Remove all meta data and veliger counts from the dataframe 
meta_drop <- c('Station_ID', 'Sample_Year', "Sample_Month","Sample_mo_yr",'ID', 
               'veliger quagga_no.L', # from zooplankton data - auto-corrected to status
               'sum_veligers_L', 'mu_veligers_L', 'sd_veligers')
wqdata_yrmo <- gbm_wqinput_yrmo[,!colnames(gbm_wqinput_yrmo) %in% meta_drop] 
wqdata_yrmo$Status <- as.factor(wqdata_yrmo$Status) # make sure status is a factor

planktondata_yrmo <- gbm_planktoninput_yrmo[,!colnames(gbm_planktoninput_yrmo) %in% meta_drop] 
planktondata_yrmo$Status <- as.factor(planktondata_yrmo$Status) # make sure status is a factor

# Adjust the column names to make them more readable and R friendly
colClean <- function(x) {colnames(x) <- gsub("mu_", "", colnames(x)) # remove mu tag
colnames(x) <- gsub("_mg.L", "", colnames(x)) # remove unit
colnames(x) <- gsub("_no.L", "", colnames(x)) # remove unit
colnames(x) <- gsub("Chla_Pheo_664.665a", "Chla/Pheophytin", colnames(x)) # change label
colnames(x) <- gsub("Deg_C", 'Temperature', colnames(x)) # replace this label
colnames(x) <- gsub("Alk", 'Alkalinity', colnames(x)) # replace this label
colnames(x) <- gsub("^Ca$", "Ca2+", colnames(x))
colnames(x) <- gsub("^K$", 'K+', colnames(x)) 
colnames(x) <- gsub("^SO4$", 'SO42-', colnames(x)) 
colnames(x) <- gsub("^F$", 'F-', colnames(x)) 
colnames(x) <- gsub(" ","_",  colnames(x))} #change spaces to underscores - more r friendly

colnames(wqdata_yrmo) <- colClean(wqdata_yrmo)
colnames(planktondata_yrmo) <- colClean(planktondata_yrmo)

# How much data are we working with?
dim(select(wqdata_yrmo,-c(Waterbody, Status))) # 179 observations, 14 WQ variables
dim(select(planktondata_yrmo,-c(Waterbody, Status))) # 139 observations, 24 variables (13 phytoplankton, 11 zooplankton)

table(gbm_wqinput_yrmo$Sample_mo_yr, gbm_wqinput_yrmo$Station_ID) #A-D_April_2021 missing
length(unique(gbm_wqinput_yrmo$ID)) # n = 179

table(gbm_planktoninput_yrmo$Sample_mo_yr, gbm_planktoninput_yrmo$Station_ID) #A-D_April_2021 missing
length(unique(gbm_planktoninput_yrmo$ID)) # n = 139

####### Section 4: Splitting the Station Year-Month data set for training and testing with caret ####

## WATER QUALTIY DATA MONTHLY MEANS *****************************************************************************
# Two forms of pre-processing must be done before splitting the data 
# or the trained model may not work with the test set
# 1. impute missing values (using the median of the whole col)
wqdata_yrmo_imp <- wqdata_yrmo
for(i in 3:ncol(wqdata_yrmo_imp)) {
  wqdata_yrmo_imp[ , i][is.na(wqdata_yrmo_imp[ , i])] <- median(wqdata_yrmo_imp[ , i], na.rm=TRUE)
}
# 2. remove columns with near zero variance - will not create splits in tree
wqdata_yrmo_nzv <- nearZeroVar(wqdata_yrmo_imp)
if(length(wqdata_yrmo_nzv > 0)) {
  wqdata_yrmo_pp <- wqdata_yrmo_imp[, -wqdata_yrmo_nzv]
}else{wqdata_yrmo_pp <- wqdata_yrmo_imp}

dim(select(wqdata_yrmo,-c(Status, Waterbody))) # row, cols before 
dim(select(wqdata_yrmo_pp,-c(Status, Waterbody))) # row, cols after

# Splitting the Bernoulli dataset into a training and testing set
# Simple Splitting the data set into a training (80%) and testing (20%) set
set.seed(2022) # constrains random sampling for reproducibility
wq_train_index_yrmo <- createDataPartition(wqdata_yrmo_pp$Status, # response var
                                           p = .8, # proportion of data to put into training set
                                           list = FALSE, 
                                           times = 1)
train_wqdata_yrmo <- wqdata_yrmo_pp[wq_train_index_yrmo, ] # put those rows in the training data set 
test_wqdata_yrmo <- wqdata_yrmo_pp[-wq_train_index_yrmo, ] # put the rest of the rows in the testing data set

# Make sure to remove waterbody from the modeling data
train_wqdata_yrmo <- train_wqdata_yrmo[,!colnames(train_wqdata_yrmo) == 'Waterbody']
test_wqdata_yrmo <- test_wqdata_yrmo[,!colnames(test_wqdata_yrmo) == 'Waterbody']

## PLANKTON DATA MONTHLY MEANS *****************************************************************************
# Two forms of pre-processing must be done before splitting the data 
# or the trained model may not work with the test set
# 1. impute missing values (using the median of the whole col)
planktondata_yrmo_imp <- planktondata_yrmo
for(i in 3:ncol(planktondata_yrmo_imp)) {
  planktondata_yrmo_imp[ , i][is.na(planktondata_yrmo_imp[ , i])] <- median(planktondata_yrmo_imp[ , i], na.rm=TRUE)
}
# 2. remove columns with near zero variance - will not create splits in tree
planktondata_yrmo_nzv <- nearZeroVar(planktondata_yrmo_imp)
if(length(planktondata_yrmo_nzv > 0)) {
  planktondata_yrmo_pp <- planktondata_yrmo_imp[, -planktondata_yrmo_nzv]
}else{planktondata_yrmo_pp <- planktondata_yrmo_imp}

dim(select(planktondata_yrmo,-c(Status, Waterbody))) # row, cols before 
dim(select(planktondata_yrmo_pp,-c(Status, Waterbody))) # row, cols after

# Splitting the Bernoulli dataset into a training and testing set
# Simple Splitting the data set into a training (80%) and testing (20%) set
set.seed(2022) # constrains random sampling for reproducibility
plankton_train_index_yrmo <- createDataPartition(planktondata_yrmo_pp$Status, # response var
                                                 p = .8, # proportion of data to put into training set
                                                 list = FALSE, 
                                                 times = 1)
train_planktondata_yrmo <- planktondata_yrmo_pp[plankton_train_index_yrmo, ] # put those rows in the training data set 
test_planktondata_yrmo <- planktondata_yrmo_pp[-plankton_train_index_yrmo, ] # put the rest of the rows in the testing data set

# Make sure to remove waterbody from the modeling data
train_planktondata_yrmo <- train_planktondata_yrmo[,!colnames(train_planktondata_yrmo) == 'Waterbody']
test_planktondata_yrmo <- test_planktondata_yrmo[,!colnames(test_planktondata_yrmo) == 'Waterbody']

#### CLASSIFICATION GRADIENT BOOSTED MACHINES (CARET) #### ----------------------------------------------------------------------------------------------
####### Section 5: Control and Tuning parameters for GBM caret ####
# https://topepo.github.io/caret/model-training-and-tuning.html
# https://www.youtube.com/watch?v=GZafoSe5DvI

# The train function can be used to
# evaluate, using resampling, the effect of model tuning parameters on performance
# choose the "optimal" model across these parameters
# estimate model performance from a training set

# We are going to be using a Gradient Boosted Model from caret which is appropriate for classification
# modelLookup('gbm') # uncomment to look at some of the details for this model such as the tuning parameters 
# note: shrinkage = learning rate

# This function gives more detailed summary of the best model
get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

# To start building our model: let's define how we want to preprocess and resample the data 
# These 'fit' objects will be argument inputs in the train() function 

# First let's define how to preprocess the data by centering, scaling, and transforming
gbmfit_preProcess <- c(#'medianImpute', # imputes missing values using median
  #'nzv', # removes near zero variance vars
  
  # medianImpute and nzv must be done before splitting the training and test sets 
  # OR the training and test sets will not line up in the predict function
  # centering and scale must be done within the model so that information 
  # does not leak from the training to testing sets
  
  'center', # subtracts mean from predictors
  'scale', # divides predictors by the standard deviation
  'YeoJohnson') # generalized Box-Cox for data with zero or negative smallest value

# Second let's define how we want to resample the data using k-fold cross validation 
gbmfit_trainControl <- trainControl(method = "repeatedcv", # aka repeated K-fold Cross-Validation
                                    number = 10, # k = 10 is common default value in machine learning
                                    repeats = 10, # repeat CV x times - suited for small data sets (usually = 3,5,10) 
                                    savePredictions = 'final', # saves predictions for optimal tuning parameter
                                    classProbs = TRUE, # should class probabilities be returned
                                    selectionFunction = 'best',# how to select tuning parameters
                                    summaryFunction=multiClassSummary) 


# Creating a tuning grid for the model
# For a gradient boosting machine (GBM) model, there are three main tuning parameters:
# number of iterations, i.e. trees, (called n.trees in the gbm function)
# complexity of the tree, called interaction.depth
# learning rate: how quickly the algorithm adapts, called shrinkage
# the minimum number of training set samples in a node to commence splitting (n.minobsinnode)

gbmGrid <- expand.grid(interaction.depth = c(1:10), # tries 1 to 10
                       n.trees = (1:6)*50, # tries 50 to 300 by jumps of 50
                       shrinkage = c(0.1,0.2,0.3), # learning rate
                       n.minobsinnode = 2) # increase this as we get more data 

####### Section 6: Station Year-Month Gradient Boosted Machine - Training and Tuning with caret; Table S6, Fig 4 ####

# MODEL EVALUATION NOTES
# Accuracy = # test water bodies correctly classified / # test water bodies
# Kappa = better 'accuracy' metric for unbalanced classes 
#         a value < 0 is indicating no agreement, 0–0.20 as slight, 0.21–0.40 as fair, 
#         0.41–0.60 as moderate, 0.61–0.80 as substantial, and 0.81–1 as almost perfect agreement

## WQ DATASET GBM MONTHLY MEANS *****************************************************************************************************************
set.seed(2022) # constrains random sampling for reproducibility
gbm_fit_wqyrmo <- caret::train(Status ~ .,
                               data = train_wqdata_yrmo, 
                               method = 'gbm', # gradient boosted machine
                               trControl = gbmfit_trainControl,
                               preProcess = gbmfit_preProcess,
                               tuneGrid = gbmGrid,
                               verbose = F) # suppresses iteration text output
# gbm_fit_wqyrmo # This will take a minute to run

# Visually inspect how these different models perform
ggplot(gbm_fit_wqyrmo)+theme_bw()

# Select the best model
# this is not simply the tree with the highest accuracy
# it is the simplest tree with high accuracy and kappa
t(get_best_result(gbm_fit_wqyrmo)) # examine that model that GBM chose
plot(varImp(gbm_fit_wqyrmo))
varImp(gbm_fit_wqyrmo)

# Extracting Predictions and Class Probabilities
predicted.raw_wqyrmo <- predict(gbm_fit_wqyrmo, newdata = test_wqdata_yrmo, type = "raw")
predicted.prob_wqyrmo <- predict(gbm_fit_wqyrmo, newdata = test_wqdata_yrmo, type = "prob")

# Compute the confusion matrix
cm_wqdata_yrmo <- caret::confusionMatrix(reference = test_wqdata_yrmo$Status, data = predicted.raw_wqyrmo, mode='everything', positive='Established');cm_wqdata_yrmo
# https://topepo.github.io/caret/measuring-performance.html#class

# Which sample is being miss-classified?
cbind(test_wqdata_yrmo[1],predicted.raw_wqyrmo) # none 

# Calculate ROC curve - two class only
roc.gbm_fit_wqyrmo <- roc(test_wqdata_yrmo$Status,predicted.prob_wqyrmo[,'Negative'])
# plot the ROC curve
plot(roc.gbm_fit_wqyrmo, col=c(3)) 

## PLANKTON DATASET GBM MONTHLY MEANS *****************************************************************************************************************
set.seed(2022) # constrains random sampling for reproducibility
gbm_fit_planktonyrmo <- caret::train(Status ~ .,
                                     data = train_planktondata_yrmo, 
                                     method = 'gbm', # gradient boosted machine
                                     trControl = gbmfit_trainControl,
                                     preProcess = gbmfit_preProcess,
                                     tuneGrid = gbmGrid,
                                     verbose = F) # suppresses iteration text output
# gbm_fit_planktonyrmo # This will take a minute to run

# Visually inspect how these different models perform
ggplot(gbm_fit_planktonyrmo)+theme_bw()

# Select the best model
# this is not simply the tree with the highest accuracy
# it is the simplest tree with high accuracy and kappa
t(get_best_result(gbm_fit_planktonyrmo)) # examine that model that GBM chose
plot(varImp(gbm_fit_planktonyrmo))
varImp(gbm_fit_planktonyrmo)

# Extracting Predictions and Class Probabilities
predicted.raw_planktonyrmo <- predict(gbm_fit_planktonyrmo, newdata = test_planktondata_yrmo, type = "raw")
predicted.prob_planktonyrmo <- predict(gbm_fit_planktonyrmo, newdata = test_planktondata_yrmo, type = "prob")

# Compute the confusion matrix
cm_planktondata_yrmo <- caret::confusionMatrix(reference = test_planktondata_yrmo$Status, data = predicted.raw_planktonyrmo, mode='everything', positive='Established'); cm_planktondata_yrmo
# https://topepo.github.io/caret/measuring-performance.html#class

# Which sample is being miss-classified?
cbind(test_planktondata_yrmo[1],predicted.raw_planktonyrmo) # B-D_July_2023

# Calculate ROC curve - two class only
roc.gbm_fit_planktonyrmo <- roc(test_planktondata_yrmo$Status,predicted.prob_planktonyrmo[,'Negative'])
# plot the ROC curve
plot(roc.gbm_fit_planktonyrmo, col=c(3)) 

# Table S6
t(get_best_result(gbm_fit_wqyrmo)) # examine that model that GBM chose
get_best_result(gbm_fit_wqyrmo)[c(1:4)]; cm_wqdata_yrmo$overall
t(get_best_result(gbm_fit_planktonyrmo)) # examine that model that GBM chose
get_best_result(gbm_fit_planktonyrmo)[c(1:4)]; cm_planktondata_yrmo$overall

# Fig 4
plot(varImp(gbm_fit_wqyrmo))
varImp(gbm_fit_wqyrmo)
plot(varImp(gbm_fit_planktonyrmo))
varImp(gbm_fit_planktonyrmo)

#### Section 7a: Figure 5a, S4a; Water Quality Partial dependence plots ####

# Get sample value and probability data for PDPs
pdp_TDS <-  partial(gbm_fit_wqyrmo, pred.var = "TDS", type = "classification", prob = T, rug = T, plot = F)
pdp_K <- partial(gbm_fit_wqyrmo, pred.var = "K+", type = "classification", prob = T, rug = T, plot = F)
pdp_F <- partial(gbm_fit_wqyrmo, pred.var = "F-", type = "classification", prob = T, rug = T, plot = F)
pdp_Total_P <- partial(gbm_fit_wqyrmo, pred.var = "TP", type = "classification", prob = T, rug = T, plot = F)
pdp_SiO2 <- partial(gbm_fit_wqyrmo, pred.var = "SiO2", type = "classification", prob = T, rug = T, plot = F)
pdp_Chla <- partial(gbm_fit_wqyrmo, pred.var = "Chla", type = "classification", prob = T, rug = T, plot = F)
pdp_SO4 <- partial(gbm_fit_wqyrmo, pred.var = "SO42-", type = "classification", prob = T, rug = T, plot = F)
pdp_Total_N <- partial(gbm_fit_wqyrmo, pred.var = "TN", type = "classification", prob = T, rug = T, plot = F)
pdp_DO <- partial(gbm_fit_wqyrmo, pred.var = "DO", type = "classification", prob = T, rug = T, plot = F)
pdp_Pheophytin <- partial(gbm_fit_wqyrmo, pred.var = "Pheophytin", type = "classification", prob = T, rug = T, plot = F)
pdp_Alk <- partial(gbm_fit_wqyrmo, pred.var = "Alkalinity", type = "classification", prob = T, rug = T, plot = F)
pdp_Deg_C <- partial(gbm_fit_wqyrmo, pred.var = "Temperature", type = "classification", prob = T, rug = T, plot = F)
pdp_Ca <- partial(gbm_fit_wqyrmo, pred.var = "Ca2+", type = "classification", prob = T, rug = T, plot = F)


# Plot Osmoregulation PDPs
tds_pdp <- ggplot() + ylim(0.45, 0.7) + theme_bw() +
  geom_path(data = pdp_TDS, aes(y = yhat, x = TDS), linewidth = 0.75, col = "#3366FF") +
  geom_point(data = train_wqdata_yrmo, aes(x = TDS), y = 0.44, shape = "|", size = 4, col = "#3366FF") +
  labs(x = "TDS (mg/L) (100)", y = "Probability of 'Established'") +  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
k_pdp <- ggplot() + ylim(0.45, 0.7) + theme_bw() +
  geom_path(data = pdp_K, aes(y = yhat, x = `K+`), linewidth = 0.75, col = "#3366FF") +
  geom_point(data = train_wqdata_yrmo, aes(x = `K+`), y = 0.44, shape = "|", size = 4, col = "#3366FF") +
  labs(x = expression(paste("K"^"+", " (mg/L) (73.1)")), y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
so4_pdp <- ggplot() + ylim(0.45, 0.7) + theme_bw() +
  geom_path(data = pdp_SO4, aes(y = yhat, x = `SO42-`), linewidth = 0.75, col = "#3366FF") +
  geom_point(data = train_wqdata_yrmo, aes(x = `SO42-`), y = 0.44, shape = "|", size = 4, col = "#3366FF") +
  labs(x = expression(paste("SO"["4"]^"2-", " (mg/L) (20.9)")), y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
grid.arrange(tds_pdp, k_pdp, so4_pdp, nrow = 1, widths = c(1, 0.8, 0.8), 
             top = grid::textGrob("Osmoregulation",  x = 0.01, hjust = 0, gp = gpar(fontsize = 14, col = "#3366FF")))

# Plot unknown relationship PDPs
f_pdp <- ggplot() + ylim(0.35, 0.7) + theme_bw() +
  geom_path(data = pdp_F, aes(y = yhat, x = `F-`), linewidth = 0.75, col = "black") +
  geom_point(data = train_wqdata_yrmo, aes(x = `F-`), y = 0.34, shape = "|", size = 4, col = "black") +
  labs(x = expression(paste("F"^"-", " (mg/L) (70.9)")), y = "Probability of 'Established'") + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), 
        plot.title = element_text(hjust = 4, size = 14))
SiO2_pdp <- ggplot() + ylim(0.35, 0.7) + theme_bw() +
  geom_path(data = pdp_SiO2, aes(y = yhat, x = SiO2), linewidth = 0.75, col = "black") +
  geom_point(data = train_wqdata_yrmo, aes(x = SiO2), y = 0.34, shape = "|", size = 4, col = "black") +
  labs(x = expression(paste("SiO"["2"], " (mg/L) (28.3)")), y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
grid.arrange(f_pdp,SiO2_pdp, nrow = 1, ncol = 2, widths = c(1, 0.8),
             top = grid::textGrob("Unknown relationships",  x = 0.01, hjust = 0,  gp = gpar(fontsize = 14)))

# Plot Phytoplankton cycling PDPs
TP_pdp <- ggplot() + ylim(0.54, 0.68) + theme_bw() +
  geom_path(data = pdp_Total_P, aes(y = yhat, x = TP), linewidth = 0.75, col = "#339966") +
  geom_point(data = train_wqdata_yrmo, aes(x = TP), y = 0.5375, shape = "|", size = 4, col = "#339966") +
  labs(x = "TP (mg/L) (41.0)", y = "Probability of 'Established'") +  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
Chla_pdp <- ggplot() + ylim(0.54, 0.68) + theme_bw() +
  geom_path(data = pdp_Chla, aes(y = yhat, x = Chla), linewidth = 0.75, col = "#339966") +
  geom_point(data = train_wqdata_yrmo, aes(x = Chla), y = 0.5375, shape = "|", size = 4, col = "#339966") +
  labs(x = "Chl a (mg/L) (28.0)", y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
grid.arrange(TP_pdp, Chla_pdp, nrow = 1, widths = c(1, 0.8), 
             top = grid::textGrob("Phytoplankton cycling",  x = 0.01, hjust = 0, gp = gpar(fontsize = 14, col = "#339966")))


# Plot <10 VarIMP PDPs for appendix
TN_pdp <- ggplot() + ylim(0.58, 0.642) + theme_bw() +
  geom_path(data = pdp_Total_N, aes(y = yhat, x = TN), linewidth = 0.75) +
  geom_point(data = train_wqdata_yrmo, aes(x = TN), y = 0.578, shape = "|", size = 4) +
  labs(x = "TN (mg/L) (9.9)", y = "Probability of 'Established'") +  
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
do_pdp <- ggplot() + ylim(0.58, 0.642) + theme_bw() +
  geom_path(data = pdp_DO, aes(y = yhat, x = DO), linewidth = 0.75) +
  geom_point(data = train_wqdata_yrmo, aes(x = DO), y = 0.578, shape = "|", size = 4) +
  labs(x = "DO (mg/L) (8.2)", y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
pheophytin_pdp <- ggplot() + ylim(0.58, 0.642) + theme_bw() +
  geom_path(data = pdp_Pheophytin, aes(y = yhat, x = Pheophytin), linewidth = 0.75) +
  geom_point(data = train_wqdata_yrmo, aes(x = Pheophytin), y = 0.578, shape = "|", size = 4) +
  labs(x = "Pheophytin (mg/L) (4.9)", y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
alk_pdp <- ggplot() + ylim(0.58, 0.642) + theme_bw() +
  geom_path(data = pdp_Alk, aes(y = yhat, x = Alkalinity), linewidth = 0.75) +
  geom_point(data = train_wqdata_yrmo, aes(x = Alkalinity), y = 0.578, shape = "|", size = 4) +
  labs(x = "Alkalinity (mg/L) (3.8)", y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
deg_c_pdp <- ggplot() + ylim(0.58, 0.642) + theme_bw() +
  geom_path(data = pdp_Deg_C, aes(y = yhat, x = Temperature), linewidth = 0.75) +
  geom_point(data = train_wqdata_yrmo, aes(x = Temperature), y = 0.578, shape = "|", size = 4) +
  labs(x = expression(paste("Temperature ("^"o","C) (0.88)")), y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
ca_pdp <- ggplot() + ylim(0.58, 0.642) + theme_bw() +
  geom_path(data = pdp_Ca, aes(y = yhat, x = `Ca2+`), linewidth = 0.75) +
  geom_point(data = train_wqdata_yrmo, aes(x = `Ca2+`), y = 0.578, shape = "|", size = 4) +
  labs(x = expression(paste("Ca"^"2+", " (mg/L) (0.84)")), y = NULL) +  
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
grid.arrange(TN_pdp, do_pdp, pheophytin_pdp, alk_pdp, deg_c_pdp, ca_pdp, nrow = 1, widths = c(1, 0.8, 0.8, 0.8, 0.8, 0.8))


#### Section 7b: Figure 5b, S4b; Plankton Partial dependence plots ####

# Get sample value and probability data for PDPs
pdp_Bosmina_longirostris <-  partial(gbm_fit_planktonyrmo, pred.var = "Bosmina_longirostris", type = "classification",  prob = T, rug = T, plot = F)
pdp_Stephanodiscus_parvus <-  partial(gbm_fit_planktonyrmo, pred.var = "Stephanodiscus_parvus", type = "classification",  prob = T, rug = T, plot = F)
pdp_nauplii <-  partial(gbm_fit_planktonyrmo, pred.var = "nauplii", type = "classification",  prob = T, rug = T, plot = F)
pdp_Scenedesmus_communis <-  partial(gbm_fit_planktonyrmo, pred.var = "Scenedesmus_communis", type = "classification",  prob = T, rug = T, plot = F)
pdp_Keratella_cochlearis <-  partial(gbm_fit_planktonyrmo, pred.var = "Keratella_cochlearis", type = "classification",  prob = T, rug = T, plot = F)
pdp_Plagioselmis_nannoplanctica <-  partial(gbm_fit_planktonyrmo, pred.var = "Plagioselmis_nannoplanctica", type = "classification",  prob = T, rug = T, plot = F)
pdp_Polyarthra_vulgaris <-  partial(gbm_fit_planktonyrmo, pred.var = "Polyarthra_vulgaris", type = "classification",  prob = T, rug = T, plot = F)
pdp_cyclopoid_copepodid <-  partial(gbm_fit_planktonyrmo, pred.var = "cyclopoid_copepodid", type = "classification",  prob = T, rug = T, plot = F)
pdp_Polyarthra_dolichoptera <-  partial(gbm_fit_planktonyrmo, pred.var = "Polyarthra_dolichoptera", type = "classification",  prob = T, rug = T, plot = F)
pdp_Keratella_americana <-  partial(gbm_fit_planktonyrmo, pred.var = "Keratella_americana", type = "classification",  prob = T, rug = T, plot = F)
pdp_Synchaeta_spp. <-  partial(gbm_fit_planktonyrmo, pred.var = "Synchaeta_spp.", type = "classification",  prob = T, rug = T, plot = F)
pdp_Chroococcus_microscopicus <-  partial(gbm_fit_planktonyrmo, pred.var = "Chroococcus_microscopicus", type = "classification",  prob = T, rug = T, plot = F)
pdp_Raphidiopsis_curvata <-  partial(gbm_fit_planktonyrmo, pred.var = "Raphidiopsis_curvata", type = "classification",  prob = T, rug = T, plot = F)
pdp_Cylindrospermopsis_raciborskii <-  partial(gbm_fit_planktonyrmo, pred.var = "Cylindrospermopsis_raciborskii", type = "classification",  prob = T, rug = T, plot = F)
pdp_Planktolyngbya_sp. <-  partial(gbm_fit_planktonyrmo, pred.var = "Planktolyngbya_sp.", type = "classification",  prob = T, rug = T, plot = F)
pdp_Chlorella_spp. <-  partial(gbm_fit_planktonyrmo, pred.var = "Chlorella_spp.", type = "classification",  prob = T, rug = T, plot = F)
pdp_Planktolyngbya_limnetica <-  partial(gbm_fit_planktonyrmo, pred.var = "Planktolyngbya_limnetica", type = "classification",  prob = T, rug = T, plot = F)
pdp_Brachionus_angularis <-  partial(gbm_fit_planktonyrmo, pred.var = "Brachionus_angularis", type = "classification",  prob = T, rug = T, plot = F)
pdp_Keratella_tropica <-  partial(gbm_fit_planktonyrmo, pred.var = "Keratella_tropica", type = "classification",  prob = T, rug = T, plot = F)
pdp_Cylindrospermopsis_sp. <-  partial(gbm_fit_planktonyrmo, pred.var = "Cylindrospermopsis_sp.", type = "classification",  prob = T, rug = T, plot = F)
pdp_Anabaenopsis_circularis <-  partial(gbm_fit_planktonyrmo, pred.var = "Anabaenopsis_circularis", type = "classification",  prob = T, rug = T, plot = F)

# Plot competitor PDPs
Bl_pdp <- ggplot() + ylim(0.35, 0.66) + theme_bw() +
  geom_path(data = pdp_Bosmina_longirostris, aes(y = yhat, x = Bosmina_longirostris), linewidth = 0.75, col = "#CC9900") +
  geom_point(data = train_planktondata_yrmo, aes(x = Bosmina_longirostris), y = 0.34, shape = "|", size = 4, col = "#CC9900") +
  labs(x = expression(paste(italic("Bosmina longirostris"), " (100)")), y = "Probability of 'Established'") +
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text = element_text(size = 12), 
        plot.title = element_text(hjust = -0.6, size = 14, color = "#CC9900"))
cc_pdp <- ggplot() + ylim(0.35, 0.66) + theme_bw() +
  geom_path(data = pdp_cyclopoid_copepodid, aes(y = yhat, x = cyclopoid_copepodid), linewidth = 0.75, col = "#CC9900") +
  geom_point(data = train_planktondata_yrmo, aes(x = cyclopoid_copepodid), y = 0.34, shape = "|", size = 4, col = "#CC9900") +
  labs(x = "Cyclopoid copepodid (26.9)", y = NULL) + 
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
grid.arrange(Bl_pdp, cc_pdp, nrow = 1, widths = c(1, 0.8), 
             top = grid::textGrob("Competitors",  x = 0.01, hjust = 0, gp = gpar(fontsize = 14, col = "#CC9900")))

# Plot community indicator for grazing PDP
Sp_pdp <- ggplot() + ylim(0.54, 0.82) + theme_bw() +
  geom_path(data = pdp_Stephanodiscus_parvus, aes(y = yhat, x = Stephanodiscus_parvus), linewidth = 0.75, col = "#FF6600") +
  geom_point(data = train_planktondata_yrmo, aes(x = Stephanodiscus_parvus), y = 0.535, shape = "|", size = 4, col = "#FF6600") +
  labs(x = expression(paste(italic("Stephanodiscus parvus"), " (95.3)")), y = "Probability of 'Established'") + 
  theme(axis.title.x = element_text(size = 12, color = "#339900"),
        axis.title.y = element_text(size = 12), 
        axis.text = element_text(size = 12))
grid.arrange(Sp_pdp, nrow = 1, ncol = 1, 
             top = grid::textGrob("Community indicator of grazing",  x = 0.01, hjust = 0, gp = gpar(fontsize = 14, col = "#FF6600")))

# Plot Prey PDPs
N_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_nauplii, aes(y = yhat, x = nauplii), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = nauplii), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = "Nauplii (75.8)", y = "Probability of 'Established'") +  
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))
Sc_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_Scenedesmus_communis, aes(y = yhat, x = Scenedesmus_communis), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = Scenedesmus_communis), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = expression(paste(italic("Scenedesmus communis"), " (72.0)")), y = NULL) +  
  theme(axis.title.x = element_text(size = 12, color = "#339900"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
Kc_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_Keratella_cochlearis, aes(y = yhat, x = Keratella_cochlearis), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = Keratella_cochlearis), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = expression(paste(italic("Keratella cochlearis")," (67.5)")), y = NULL) +  
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
Pn_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_Plagioselmis_nannoplanctica, aes(y = yhat, x = Plagioselmis_nannoplanctica), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = Plagioselmis_nannoplanctica), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = expression(paste(italic("Plagioselmis nannoplanctica"), " (43.8)")), y = NULL) +  
  theme(axis.title.x = element_text(size = 12, color = "#339900"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
Pv_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_Polyarthra_vulgaris, aes(y = yhat, x = Polyarthra_vulgaris), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = Polyarthra_vulgaris), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = expression(paste(italic("Polyarthra vulgaris"), " (38.7)")), y = "Probability of 'Established'") +  
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))
Pd_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_Polyarthra_dolichoptera, aes(y = yhat, x = Polyarthra_dolichoptera), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = Polyarthra_dolichoptera), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = expression(paste(italic("Polyarthra dolichoptera")," (24.7)")), y = NULL) +  
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
Ka_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_Keratella_americana, aes(y = yhat, x = Keratella_americana), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = Keratella_americana), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = expression(paste(italic("Keratella americana"), " (20.1)")), y = NULL) +  
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
S_pdp <- ggplot() + ylim(0.50, 0.82) + theme_bw() +
  geom_path(data = pdp_Synchaeta_spp., aes(y = yhat, x = Synchaeta_spp.), linewidth = 0.75, col = "#CC0066") +
  geom_point(data = train_planktondata_yrmo, aes(x = Synchaeta_spp.), y = 0.495, shape = "|", size = 4, col = "#CC0066") +
  labs(x = expression(paste(italic("Synchaeta spp."), " (15.0)")), y = NULL) +  
  theme(axis.title = element_text(size = 12, color = "#663300"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_blank())
grid.arrange(N_pdp, Sc_pdp, Kc_pdp, Pn_pdp, 
             Pv_pdp, Pd_pdp, Ka_pdp, S_pdp, nrow = 2, ncol = 4, widths = c(1, 0.8, 0.8, 0.8),
             top = grid::textGrob("Prey",  x = 0.01, hjust = 0, gp = gpar(fontsize = 14, col = "#CC0066")))

# Plot <10 VarIMP PDPs for appendix
Cm_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Chroococcus_microscopicus, aes(y = yhat, x = Chroococcus_microscopicus), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Chroococcus_microscopicus), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Chroococcus microscopicus"), " (9.6)")), y = "Probability of 'Established'") +  
  theme(axis.title.x = element_text(size = 10, color = "#339900"),
        axis.title.y = element_text(size = 10), 
        axis.text = element_text(size = 10))
Rc_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Raphidiopsis_curvata, aes(y = yhat, x = Raphidiopsis_curvata), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Raphidiopsis_curvata), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Raphidiopsis curvata")," (8.8)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#339900"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank())
Cr_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Cylindrospermopsis_raciborskii, aes(y = yhat, x = Cylindrospermopsis_raciborskii), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Cylindrospermopsis_raciborskii), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Cylindrospermopsis raciborskii")," (7.2)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#339900"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank())
P_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Planktolyngbya_sp., aes(y = yhat, x = Planktolyngbya_sp.), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Planktolyngbya_sp.), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Planktolyngbya sp.")," (5.8)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#339900"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank())
C_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Chlorella_spp., aes(y = yhat, x = Chlorella_spp.), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Chlorella_spp.), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Chlorella spp.")," (3.9)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#339900"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank())
Pl_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Planktolyngbya_limnetica, aes(y = yhat, x = Planktolyngbya_limnetica), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Planktolyngbya_limnetica), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Planktolyngbya limnetica"), " (3.2)")), y = "Probability of 'Established'") +  
  theme(axis.title.x = element_text(size = 10, color = "#339900"),
        axis.title.y = element_text(size = 10), 
        axis.text = element_text(size = 10))
Ba_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Brachionus_angularis, aes(y = yhat, x = Brachionus_angularis), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Brachionus_angularis), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Brachionus angularis")," (2.9)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#663300"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank())
Kt_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Keratella_tropica, aes(y = yhat, x = Keratella_tropica), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Keratella_tropica), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Keratella tropica")," (0.92)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#663300"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank())
Cy_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Cylindrospermopsis_sp., aes(y = yhat, x = Cylindrospermopsis_sp.), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Cylindrospermopsis_sp.), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Cylindrospermopsis sp.")," (0.77)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#339900"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank())
Ac_pdp <- ggplot() + ylim(0.597, 0.599) + theme_bw() +
  geom_path(data = pdp_Anabaenopsis_circularis, aes(y = yhat, x = Anabaenopsis_circularis), linewidth = 0.75) +
  geom_point(data = train_planktondata_yrmo, aes(x = Anabaenopsis_circularis), y = 0.59695, shape = "|", size = 4) +
  labs(x = expression(paste(italic("Anabaenopsis circularis")," (0.18)")), y = NULL) +  
  theme(axis.title = element_text(size = 10, color = "#339900"),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_blank()) 
grid.arrange(Cm_pdp, Rc_pdp, Cr_pdp, P_pdp, C_pdp, Pl_pdp, Ba_pdp, Kt_pdp, Cy_pdp, Ac_pdp, nrow = 2, ncol = 5, widths = c(1, 0.8, 0.8, 0.8, 0.8))


#### Section 8: Figure 6, 7, S5, S6; Additional Post hoc variable investigations ####

# Some variables appear to have high value "outliers" - make sure these are not all from the same station/waterbody/date
all_pp_vars <- merge(wqdata_yrmo_pp, planktondata_yrmo_pp, by = c("row.names", "Status", "Waterbody"))
rownames(all_pp_vars) <- all_pp_vars$Row.names
all_pp_vars$Row.names <- NULL
all_pp_vars <- all_pp_vars[,-c(1:2)]
max_values <- data.frame(Var = character(), Sample_ID = character(), stringsAsFactors = FALSE)

# Iterate over each column
for (col in colnames(all_pp_vars)) {  # Exclude the first column (Sample_ID)
  max_val <- max(all_pp_vars[, col])  # Find the maximum value in the column
  sample_id <- rownames(all_pp_vars)[which(all_pp_vars[, col] == max_val)]  # Find the row name corresponding to the maximum value
  max_values <- rbind(max_values, data.frame(Var = col, Sample_ID = sample_id))  # Append to the results dataframe
}
table(max_values$Sample_ID) # no clear outlier sample(s) that would drive patterns in pdp
max_values[max_values$Sample_ID == "S-BC_April_2021",] # except this sample id has a bunch of max values

## Fig. S5
# Raw WQ sample density plots for all waterbodies
wqdata_yrmo %>% 
 pivot_longer(-c(Waterbody, Status)) %>% 
 mutate(Waterbody = as.factor(Waterbody), Status = as.factor(Status)) %>% 
 ggplot(aes(x = value)) + 
 geom_density(aes(col = Waterbody, linetype = Status), linewidth = 0.75) +
 geom_point(aes(col = Waterbody), y = 0, shape = "|", size = 3) +
 facet_wrap(~name, scales = 'free') + theme_bw() +
 labs(col = "Waterbody",linetype = "Status") +  
 theme(axis.title.y.left = element_text(size = 10),
       axis.title.y.right = element_text(size = 10),
       axis.title.x = element_text(size = 10), 
       strip.text = element_text(size=10))

## Fig. S6
# Raw plankton sample density plots for all waterbodies (ln(x+1)) 
planktondata_yrmo %>% 
 pivot_longer(-c(Waterbody, Status)) %>% 
 mutate(Waterbody = as.factor(Waterbody), Status = as.factor(Status)) %>% 
 ggplot(aes(log(value+1))) + 
 geom_density(aes(col = Waterbody, linetype = Status), linewidth = 0.75) +
 geom_point(aes(col = Waterbody),y = 0, shape = "|", size = 3) +
 facet_wrap(~name, scales = 'free') + theme_bw() +
 labs(col = "Waterbody",linetype = "Status") +  
 theme(axis.title.y.left = element_text(size = 10),
       axis.title.y.right = element_text(size = 10),
       axis.title.x = element_text(size = 10), 
       strip.text = element_text(size=10))

 
# Plankton sample density plots - Established waterbodies pooled together
wq_sample_density <- wqdata_yrmo
wq_sample_density$New_waterbody <- ifelse(wq_sample_density$Status == "Established", "Established waterbodies", wq_sample_density$Waterbody)
wq_sample_density$New_waterbody <- factor(wq_sample_density$New_waterbody,
                                                   levels = c("Established waterbodies",
                                                              "Bartlett Reservoir", "Theodore Roosevelt Lake"))

## Fig. 6
wq_sample_density %>% 
  pivot_longer(-c(Waterbody, New_waterbody, Status)) %>% 
  mutate(New_waterbody = as.factor(New_waterbody)) %>% 
  ggplot(aes(x = value)) + 
  geom_density(aes(col = New_waterbody), linewidth = 0.75) +
  geom_point(aes(col = New_waterbody), y = 0, shape = "|", size = 3) +
  scale_color_manual(values = c('lightcoral','aquamarine4','mediumturquoise')) +
  facet_wrap(~name, scales = 'free') + theme_bw() +
  labs(col = "Waterbody") +  
  theme(axis.title.y.left = element_text(size = 10),
        axis.title.y.right = element_text(size = 10),
        axis.title.x = element_text(size = 10), 
        strip.text = element_text(size=10))

 
# Plankton sample density plots - Established waterbodies pooled together  
plankton_sample_density <- planktondata_yrmo 
plankton_sample_density$New_waterbody <- ifelse(plankton_sample_density$Status == "Established", "Established waterbodies", plankton_sample_density$Waterbody)
plankton_sample_density$New_waterbody <- factor(plankton_sample_density$New_waterbody, 
                                    levels = c("Established waterbodies", 
                                               "Bartlett Reservoir", "Theodore Roosevelt Lake"))
## Fig. 7
plankton_sample_density %>% 
  pivot_longer(-c(Waterbody, New_waterbody, Status)) %>% 
  mutate(New_waterbody = as.factor(New_waterbody)) %>% 
  ggplot(aes(x = log(value+1))) + 
  geom_density(aes(col = New_waterbody), linewidth = 0.75) +
  geom_point(aes(col = New_waterbody), y = 0, shape = "|", size = 3) +
  scale_color_manual(values = c('lightcoral','aquamarine4','mediumturquoise')) +
  facet_wrap(~name, scales = 'free') + theme_bw() +
  labs(col = "Waterbody") +  
  theme(axis.title.y.left = element_text(size = 10),
        axis.title.y.right = element_text(size = 10),
        axis.title.x = element_text(size = 10), 
        strip.text = element_text(size=10))

# https://stats.stackexchange.com/questions/215270/my-density-plot-in-r-has-values-beyond-1-how-can-i-fix-this
