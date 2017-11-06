#! /usr/bin/env Rscript
rm(list = ls(all = TRUE)) #CLEAR WORKSPACE
#Load Required Packages
library('caret')
library('glmnet')

############################
# Load the Data, choose target, create train and test sets
############################

Data <- read.csv("/home/binyang_ni/Feature_Engineering/Kaggle_method/overfitting.csv", header=TRUE)

#Choose Target
Data$Target <- as.factor(ifelse(Data$Target_Practice ==1,'X1','X0'))
Data$Target_Evaluate = NULL
Data$Target_Leaderboard = NULL
Data$Target_Practice = NULL
xnames <- setdiff(names(Data),c('Target','case_id','train'))

#Order
Data <- Data[,c('Target','case_id','train',xnames)]

#Split to train and test
trainset = Data[Data$train == 0,]
testset = Data[Data$train == 1,]

#Remove unwanted columns
trainset$case_id = NULL
trainset$train = NULL


####################################
# RFE parameters
####################################
library(ipred)
library(e1071)

#Custom Functions
glmnetFuncs <- caretFuncs #Default caret functions

glmnetFuncs$summary <-  twoClassSummary

glmnetFuncs$rank <- function (object, x, y) {
  vimp <- sort(object$finalModel$beta[, 1])
  vimp <- as.data.frame(vimp)
  vimp$var <- row.names(vimp)
  vimp$'Overall' <- seq(nrow(vimp),1)
  vimp
}

MyRFEcontrol <- rfeControl(
  functions = glmnetFuncs,
  method = "boot",
  number = 25,
  rerank = FALSE,
  returnResamp = "final",
  saveDetails = FALSE,
  verbose = TRUE)


####################################
# Training parameters
####################################
MyTrainControl=trainControl(
  method = "boot",
  number=25,
  returnResamp = "all",
  classProbs = TRUE,
  summaryFunction=twoClassSummary
)

####################################
# Setup Multicore
####################################
#source:
#http://www.r-bloggers.com/feature-selection-using-the-caret-package/
if ( require("multicore", quietly = TRUE, warn.conflicts = FALSE) ) {
  MyRFEcontrol$workers <- multicore:::detectCores()
  MyRFEcontrol$computeFunction <- mclapply
  MyRFEcontrol$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
  
  MyTrainControl$workers <- multicore:::detectCores()
  MyTrainControl$computeFunction <- mclapply
  MyTrainControl$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
}


####################################
# Select Features-GLMNET
####################################

x <- trainset[,xnames]
y <- trainset$Target

RFE <- rfe(x,y,sizes = seq(50,200,by=10),
           metric = "ROC",maximize=TRUE,rfeControl = MyRFEcontrol,
           method='glmnet',
           tuneGrid = expand.grid(.alpha=0,.lambda=c(0.01,0.02)),
           trControl = MyTrainControl)

NewVars <- RFE$optVariables
RFE
plot(RFE)

FL <- as.formula(paste("Target ~ ", paste(NewVars, collapse= "+"))) #RFE






####################################
# Fit a GLMNET Model
####################################

model <- train(FL,data=trainset,method='glmnet',
               metric = "ROC",
               tuneGrid = expand.grid(.alpha=c(0,1),.lambda=seq(0,.25,by=0.005)),
               trControl=MyTrainControl)
model
plot(model, metric='ROC')
test <- predict(model, newdata=testset, type  = "prob")
colAUC(test, testset$Target)

predictions <- test

########################################
#Generate a file for submission
########################################
testID  <- testset$case_id
submit_file = data.frame('Zach'=predictions[,1])
write.csv(submit_file, file="AUC_ZACH.txt", row.names = FALSE)
















