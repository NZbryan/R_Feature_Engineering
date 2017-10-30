library(glmnet)
library(caret)
library(data.table)
library(stringr)
#
## input data ##

data_pre <- function(data_input,label_input){
  # input_dat = fread('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt',stringsAsFactors=FALSE, check.names=TRUE ,sep = '\t')
  # input_dat = read.csv('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt',header = T,
  #                      row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t')
  # lable_info = fread('/home/binyang_ni/R/R_data/glm_project_data/labe_info.txt',stringsAsFactors=FALSE, check.names=TRUE)
  input_dat = read.csv(data_input,header = T,
                       row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t')
  lable_info = fread(label_input,stringsAsFactors=FALSE, check.names=TRUE)
  
  colnames(lable_info) = c('SampleID','IsTraining','label')
  input_dat = as.data.frame(t(input_dat))
  # lable_info$label[lable_info$label=="Malignant"] = 1
  # lable_info$label[lable_info$label=="Benign"] = 0
  # lable_info$label = as.integer(lable_info$label)
  
  lable_info_train = subset(lable_info, lable_info$IsTraining==1)
  lable_info_test = subset(lable_info, lable_info$IsTraining==0)
  
  train_set = input_dat[lable_info_train$SampleID,]
  test_set = input_dat[lable_info_test$SampleID,]
  
  train_set$label = lable_info_train$label
  test_set$label = lable_info_test$label
  return(list(train_set=train_set,test_set=test_set ))
}
  

# model select par
select_par <- function(data_input,label_input, set_fold = 3, n_times = 50,my_step = 0.01){
  
  if(nchar(dirname(data_input))==1){
    my_path = getwd()
  }else{
    my_path = dirname(data_input) # get the file path
  }
  
  data_preoutput_output = data_pre(data_input,label_input) #  call the function : data_pre
  df = data_preoutput_output$train_set # return the values from data_pre
  set.seed(20170324)
  cv_Folds = createFolds(1:dim(df)[1],k=set_fold) # set the validation folds
  
  k=3
  train_dat = df[-cv_Folds[[k]],]
  vali_dat = df[cv_Folds[[k]],]
  
  train_x = as.matrix(subset(train_dat,select = -label)) # train set
  vali_x = as.matrix(subset(vali_dat,select = -label)) # validation set
  
  train_y = train_dat$label
  vali_y = vali_dat$label
  
  cv_value_error_50_count = 0
  cv_value_lambda_50_count = 0
  cv_value_df_50_count = 0
  # n_times = 50
  # i_set = seq(0,1,0.01)
  i_set = seq(0,1,my_step)
  for(i in 1:n_times){
    cv_alpha = c()
    cv_lambda_min = c()
    cv_value_error = c()
    cv_df = c()
    for(i in i_set){
      fit_cv4 = cv.glmnet(train_x, train_y, family="binomial",alpha=i,type.measure = "auc") # run the model
      cv_alpha = append(cv_alpha,i) # get the value of alpha
      cv_lambda_min = append(cv_lambda_min,fit_cv4$lambda.min) # get the min lambda
      cv_value_error = append(cv_value_error,min(fit_cv4$cvm)) # Statistical error
      cv_df = append(cv_df,fit_cv4$glmnet.fit$df[which(fit_cv4$glmnet.fit$lambda==fit_cv4$lambda.min)])# Statistical df
    }
    
    find_parameter = data.frame(cv_alpha,cv_lambda_min,cv_value_error,cv_df)
    cv_value_error_50_count = cv_value_error_50_count + find_parameter$cv_value_error # sum the error
    cv_value_lambda_50_count = cv_value_lambda_50_count + find_parameter$cv_lambda_min # sum the lambda
    cv_value_df_50_count = cv_value_df_50_count + find_parameter$cv_df # sum the df
  } 
  
  cv_value_error_50_mean = cv_value_error_50_count/n_times # calculatethe the average of error by the n_times
  cv_value_lambda_50_mean = cv_value_lambda_50_count/n_times # calculatethe the average of lambda by the n_times
  cv_value_df_50_mean = cv_value_df_50_count/n_times # calculatethe the average of df by the n_times
  
  output = data.frame(cv_alpha,cv_value_error_50_count,cv_value_error_50_mean,cv_value_lambda_50_count,cv_value_lambda_50_mean,cv_value_df_50_mean)
  colnames(output) = c('alpha','error_accumulation','error_mean','lambda_accumulation','lambda_mean','df_mean') # rename  the output data.frame
  select_pra = subset(output, error_mean == min(error_mean))
  
  best_alpha = select_pra$alpha
  best_lambda = select_pra$lambda_mean
  # return(output)
  return(list(best_alpha=best_alpha, best_lambda=best_lambda, select_pra = select_pra,output=output)) #  return alpha、 lambda 、 and other relative values
}


### 

# outp3 = select_par('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt','/home/binyang_ni/R/R_data/glm_project_data/labe_info.txt',
#                   n_times = 2,my_step = 0.05)


