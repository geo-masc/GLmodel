# --------------------------------------------------------------------------- #

# Random Forest Regression for the modelling of (grassland) parameters based on spectral input data.
# An adjusted version of this script was used for the analysis presented in the PFG contribution: 
# Schwieder et al. (in review): Estimating grassland parameters from Sentinel-2: A model comparison study.
#
#   Author:  Marcel Schwieder
#   Email:   marcel.schwieder@geo.hu-berlin.de
#   Date:    02/06/2020
#
# --------------------------------------------------------------------------- #

# load required packages

library(randomForest)
library(raster)
library(ggplot2)
library(gridExtra)

# Define function to split data in training and validation
# adapted from https://gist.github.com/stephenturner/834760
# documentation: https://gettinggeneticsdone.blogspot.com/2011/02/split-data-frame-into-testing-and.html
  
splitdf <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/n))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}


### Data import

# Import field data intended for modelling as data frame.
# Tables should follow the following format:
# Rows: Samples (one field observation per line)
# Columns: Spectral bands + response Variable

InputData = read.table("*.csv")

# load raster image and crop to desired extent

RS_rst <- stack("*.TIF")

# crop image to study area
UL <- # upper left coordinate
UR <- # upper right coordinate
LL <- # lower left coordinate
LR <- # lower right coordinate
cropmask <- extent(c(UL, UR, LL, LR))
RS_rst_cropped <- crop(RS_rst, cropmask)

# Convert raster to data frame and re-scale reflectance values if necessary
rstDF <- as.data.frame(RS_rst_cropped)
rstDF <- rstDF/10000

# change the raster band names
clNames <- colnames(InputData[1:9]) # i.e. band names of df
colnames(rstDF) <- clNames


### Initialize data frame to store data, models and evaluation metrics per iteration

# define split of data into training and validation
n <-(10/3)

# empty matrix for storing model performance
RFperformance <- matrix(nrow = 100, ncol = 6)

# statistical measured of the regression
colnames(RFperformance) <- c("R2", "RMSE", "relRMSE_mean", "relRMSE_range", "relRMSE_quant", "expVar") 

# empty list for raster prediction
RSTpredictions <- list()

# empty list for storing RF model training and validation data
trainData <- list()
validData <- list()
RFmodel <- list()


### Loop of Random Forest Regression 

# If the response variable is in the last column varPred = 0, else numerical value that indicates the position of the variable from the last column 
# (e.g. varPred = 1 if the response variable is in the second last column)
varPred <- 0

# Define the last columns that should be considered as prediction variables as seen from the end of the table (e.g., subVar = 1 if all but the last 
# column should be used for prediction)
subVar <- 1


for (i in 1:100){
    
    # print current status of the loop and RF model
    cat("Starting iteration ", i, "using",  tail(names(InputData), 1), "/n")
    
    # split data in training and validation data 
    splitData <- splitdf(InputData)
    valid <- splitData$trainset
    train <- splitData$testset
    
    # tune Random Forest Regression  
    rfrModel <- tuneRF(train[,1:(ncol(train)-subVar)], (train[,ncol(train)-varPred]), ntreeTry =500, stepFactor=1.5, improve=0.02, trace=TRUE, doBest=TRUE, plot=F, importance=TRUE)
    
    # predict RF model
    predict.val.rf <- predict(rfrModel, valid[,1:(ncol(valid)-subVar)])
    
    # evaluate RF model 
    r_square.rf <-  cor(valid[,ncol(valid)-varPred], predict.val.rf)^2
    RMSE.rf  <-  sqrt((sum((predict.val.rf-valid[,ncol(valid)-varPred])^2))/length(predict.val.rf))
    varExp <- rfrModel$rsq[500]*100
    rfrModel$importance

    # save model performance 
    RFperformance[i,1] <- r_square.rf
    RFperformance[i,2] <- RMSE.rf
    RFperformance[i,3] <- RMSE.rf/mean(valid[,ncol(valid)-varPred])
    RFperformance[i,4] <- RMSE.rf/ diff(range(valid[,ncol(valid)-varPred]))
    RFperformance[i,5] <- RMSE.rf/(quantile(valid[,ncol(valid)-varPred], 0.75)-quantile(valid[,ncol(valid)-varPred], 0.25))
    RFperformance[i,6] <- varExp

    # save model parameter 
    trainData[[i]] <- train
    validData[[i]] <- valid
    RFmodel[[i]] <- rfrModel
    
    # plot RF Model
    plot(rfrModel)
    
    # apply RF model to raster df
    predict_RST <- predict(rfrModel, rstDF)
    # create raster variable based on one image band 
    RSTpred <- RS_rst_cropped[[1]]
    # assign values to raster
    values(RSTpred) <- predict_RST
    # save raster prediction
    RSTpredictions[[i]] <- RSTpred
    
}


### Create mean and std RF prediction raster

# stack all raster predictions
predStack <- stack(RSTpredictions)

# compute mean raster prediction 
finalRST <- raster::calc(predStack, fun = mean)

# compute standard deviation (sd) of all raster prediction
finalRST_std <- raster::calc(predStack, fun = sd)

# save mean raster
writeRaster(finalRST, filename=paste0("/Raster/", Sensor, "_pred_mean_", varPredName, ".tif"), format="GTiff", overwrite=TRUE)

# save sd (std) raster
writeRaster(finalRST_std, filename=paste0("/Raster/", Sensor, "_pred_std_", varPredName, ".tif"), format="GTiff", overwrite=TRUE)



### end
