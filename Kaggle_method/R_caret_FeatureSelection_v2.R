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


train_x <- as.matrix(trainset[,xnames])
train_y <- trainset$Target
####################################
# test lasso
####################################

# fit_cv = cv.glmnet(train_x_sample, train_y_sample, family="binomial",alpha=1,nfolds=10, type.measure = "auc") # run the model and get the parameter lambda.1se
# fit_model = glmnet(train_x_sample, train_y_sample, family="binomial",
#                    alpha=1, lambda =fit_cv$lambda.1se) # model
# 
# feature_important = as.data.frame(as.matrix(fit_model$beta))
# feature_important = feature_important[colnames(train_x_sample),]


library(sampling)
library(methods)

feature_selection <- function(train_x,train_y,sampling_times = 50){
  beta_cal = data.frame(feature_name=colnames(train_x))
  train_y_sort = data.frame(rank=seq(length(train_y)),label=train_y)
  train_y_sort = train_y_sort[order(train_y),]
  train_x = train_x[train_y_sort$rank,]
  train_y = train_y_sort$label
  for (i in seq_len(sampling_times)){
    set_myseed = as.integer(paste(1017,i,sep = ""))
    set.seed(set_myseed)
    # sampling
    #sub1 =  sample(nrow(train_x), as.integer(nrow(train_x)*0.75), replace = F)
    sub1 = strata(train_y_sort, stratanames="label",
                  size = as.integer(c(table(train_y)[1]*0.75,table(train_y)[2]*0.75)), method =  "srswor")
    train_x_sample = train_x[sub1$ID_unit,]
    train_y_sample = train_y[sub1$ID_unit]
    
    # tuning parameter by 10 fold cv
    # glmnet(train_x_sample, train_y_sample, family="binomial",alpha=1);stop()
    fit_cv = cv.glmnet(train_x_sample, train_y_sample, family="binomial",alpha=1,nfolds=10, type.measure = "deviance") # run the model and get the parameter lambda.1se
    
    tLL <- fit_cv$glmnet.fit$nulldev - deviance(fit_cv$glmnet.fit)
    k <- fit_cv$glmnet.fit$df
    n <- fit_cv$glmnet.fit$nobs
    AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
    
    getmin = function(lambda,AIC){
      AICmin=min(AIC,na.rm=TRUE)
      idmin=AIC<=AICmin
      lambda_selected=max(lambda[idmin],na.rm=TRUE)
      lambda_selected
    }
    
    # lambda_selected = getmin(fit_cv$glmnet.fit$lambda[2:100],AICc[2:100])
    lambda_selected = getmin(fit_cv$glmnet.fit$lambda,AICc)
    # fit model
    fit_model = glmnet(train_x_sample, train_y_sample, family="binomial",
                       alpha=1, lambda =lambda_selected) # model
    
    feature_important = as.data.frame(as.matrix(fit_model$beta))
    feature_important = feature_important[colnames(train_x),]
    beta_cal = cbind(beta_cal,feature_important)
    colnames(beta_cal)[i+1] = paste('time_',i,sep = "")
    
  }
  
  return(beta_cal)
}

#output = feature_selection(data.matrix(train_x),train_y,sampling_times = 500)
output = feature_selection(train_x,train_y,sampling_times = 50)
# write.table(output,file=args[3],sep="\t",quote=F)
count_feature = data.frame(feature_name=output$feature_name,count=apply(output[2:ncol(output)]!=0,1,sum))



