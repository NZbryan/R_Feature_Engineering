#! /usr/bin/env Rscript
rm(list = ls(all = TRUE)) #CLEAR WORKSPACE
#Load Required Packages
library('caret')
library('glmnet')

############################
# Load the Data, choose target, create train and test sets
############################
x <- trainset
y <- trainset$Target

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


RFE <- rfe(x,y,sizes = seq(10,50,by=5),
           metric = "ROC",maximize=TRUE,rfeControl = MyRFEcontrol,
           method='glmnet',
           tuneGrid = expand.grid(.alpha=0,.lambda=c(0.001,0.002)),
           trControl = MyTrainControl)

NewVars <- RFE$optVariables
RFE
plot(RFE)

FL <- as.formula(paste("Target ~ ", paste(NewVars, collapse= "+"))) #RFE


