library(glmnet)
library(caret)
library(data.table)
library(pROC)

# Data preprocessing： return training set、testing set----------------------------------------------------------------
data_pre <- function(data_input,label_input){
  # input_dat = fread('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt',stringsAsFactors=FALSE, check.names=TRUE ,sep = '\t')
  # input_dat = read.csv('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt',header = T,
  #                      row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t')
  # lable_info = fread('/home/binyang_ni/R/R_data/glm_project_data/labe_info.txt',stringsAsFactors=FALSE, check.names=TRUE)
  # input_dat = read.csv(data_input,header = T,
  #                      row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t')
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

# model select parameters：return alpha、 lambda 、 and other relative values -----------------------------------------
select_par <- function(data_input,label_input, set_fold = 3, n_times = 50,my_step = 0.01){
  
  if(nchar(dirname(data_input))==1){
    my_path = getwd()
  }else{
    my_path = dirname(data_input) # get the file path
  }
  
  data_preoutput_output = data_pre(data_input,label_input) #  call the function : data_pre
  df = data_preoutput_output$train_set # return the values from data_pre
  set.seed(20170325)
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
      fit_cv4 = cv.glmnet(train_x, train_y, family="binomial",
                          alpha=i,type.measure = "auc") # run the model
      cv_alpha = append(cv_alpha,i) # get the value of alpha
      cv_lambda_min = append(cv_lambda_min,fit_cv4$lambda.min) # get the min lambda
      cv_value_error = append(cv_value_error,min(fit_cv4$cvm)) # Statistical error
      cv_df = append(cv_df,
                     fit_cv4$glmnet.fit$df[which(fit_cv4$glmnet.fit$lambda==fit_cv4$lambda.min)])# Statistical df
    }
    
    find_parameter = data.frame(cv_alpha,cv_lambda_min,cv_value_error,cv_df)
    cv_value_error_50_count = cv_value_error_50_count + find_parameter$cv_value_error # sum the error
    cv_value_lambda_50_count = cv_value_lambda_50_count + find_parameter$cv_lambda_min # sum the lambda
    cv_value_df_50_count = cv_value_df_50_count + find_parameter$cv_df # sum the df
  } 
  
  cv_value_error_50_mean = cv_value_error_50_count/n_times # calculatethe the average of error by the n_times
  cv_value_lambda_50_mean = cv_value_lambda_50_count/n_times # calculatethe the average of lambda by the n_times
  cv_value_df_50_mean = cv_value_df_50_count/n_times # calculatethe the average of df by the n_times
  
  output = data.frame(cv_alpha,cv_value_error_50_count,cv_value_error_50_mean,
                      cv_value_lambda_50_count,cv_value_lambda_50_mean,cv_value_df_50_mean)
  colnames(output) = c('alpha','error_accumulation','error_mean',
                       'lambda_accumulation','lambda_mean','df_mean') # rename  the output data.frame
  select_pra = subset(output, error_mean == min(error_mean))
  
  best_alpha = select_pra$alpha
  best_lambda = select_pra$lambda_mean
  # return(output)
  return(list(best_alpha=best_alpha, best_lambda=best_lambda, select_pra = select_pra,
              output=output)) #  return alpha、 lambda 、 and other relative values
}

# output the cross validation ROC curve image -------------------------------------------------------------------------
CV_threshold <- function(data_input,label_input,best_alpha = 0.7,best_lambda=0.01){
  if(nchar(dirname(data_input))==1){
    my_path = getwd()
  }else{
    my_path = dirname(data_input)
  }
  #
  data_preoutput_output = data_pre(data_input,label_input) #data input
  train_set = data_preoutput_output$train_set
  
  pdf_path = paste0(my_path,'/auc_plot.pdf') # pdf path
  pdf(pdf_path,width =8.5 , height = 20)
  all_thresholds = c()
  auc_threshold_list =list()
  my_step = 9
  for(i in 0:my_step){
    set_myseed = as.integer(paste(20170323,i,sep = ""))
    set.seed(set_myseed)
    cv_Folds = createFolds(1:dim(train_set)[1],k=3)
    par(mfrow = c(3,1))
    
    for(k in 1:length(cv_Folds))
    {
      train_dat = train_set[-cv_Folds[[k]],]
      vali_dat = train_set[cv_Folds[[k]],]
      
      train_x = as.matrix(subset(train_dat,select = -label))
      vali_x = as.matrix(subset(vali_dat,select = -label))
      
      train_y = train_dat$label
      vali_y = vali_dat$label
      
      # fit_cv4 = cv.glmnet(train_x, train_y, family="binomial",alpha=0.3)
      fit_bi = glmnet(train_x, train_y, family="binomial",
                      alpha=best_alpha, lambda = best_lambda) # model
      pred_vali= as.data.frame(predict(fit_bi,
                                       vali_x,type ='response'))$s0 #   predict validation data
      
      my_level = unique(train_y)
      vali_th = roc(vali_y, pred_vali,levels = my_level,ci=TRUE, 
                    of="thresholds",thresholds="best") # best thresholds
      best_thresholds = round(as.numeric(rownames(vali_th$ci$specificity)),3)
      AUC = vali_th$auc[1]
      vali_thr_auc = data.frame(best_thresholds,AUC)
      my_underline = paste0('-------------------------------',k-1,i,'-----------------------------------')
      all_thresholds = append(all_thresholds,best_thresholds)
      
      vali_thr_se_sp = data.frame(vali_th$thresholds,
                                  vali_th$sensitivities,vali_th$specificities) # sp、 se 、 threshold 
      colnames(vali_thr_se_sp) = c('thresholds' ,'sensitivities ','specificities')
      mergeouy_output_dat = list(my_underline,vali_thr_auc,vali_thr_se_sp)
      auc_threshold_list = append(auc_threshold_list,mergeouy_output_dat) # list the output
      txt_path = paste0(my_path,'/auc_thresholds.txt')
      # capture.output(mergeouy_output_dat, file = txt_path,append = TRUE)
      library(pROC)
      main_my = paste("Confidence interval of a threshold: ",'fold',k,paste('----',i+1,'th',sep=''))
      par(cex=1.2)
      plot.roc(vali_y, pred_vali,
               
               main=main_my, percent=TRUE,
               
               ci=TRUE, of="thresholds", # compute AUC (of threshold)
               
               thresholds="best", # select the (best) threshold
               
               print.thres="best",# also highlight this threshold on the plot
               
               col="#008600",
               print.auc=TRUE) 
    }
  }
  dev.off()
  
  mean_thresholds = mean(all_thresholds)
  times_count = (my_step+1)*length(cv_Folds)
  sprintf('mean of thresholds(%i times): %f' ,times_count,mean_thresholds)
  
}

# train the model with the best parameters: output a txt、a image and a model file-------------------------------------
train_model <- function(data_input,label_input,best_alpha = 0.7,best_lambda=0.01){
  if(nchar(dirname(data_input))==1){
    my_path = getwd()
  }else{
    my_path = dirname(data_input)
  }
  #
  data_preoutput_output = data_pre(data_input,label_input) #  call the function : data_pre
  train_dat = data_preoutput_output$train_set # return the values from data_pre
  
  train_x = as.matrix(subset(train_dat,select = -label)) # train set
  train_y = train_dat$label
  fit_model = glmnet(train_x, train_y, family="binomial",alpha=best_alpha, lambda = best_lambda) # model
  # save model
  my_model_path = paste0(my_path,"/glmnet_model.rds")
  saveRDS(fit_model, file=my_model_path)
  pred_train_x= as.data.frame(predict(fit_model,train_x,type ='response'))$s0
  # save txt
  my_level = unique(train_y)
  train_th = roc(train_y, pred_train_x,levels = my_level,ci=TRUE, 
                 of="thresholds",thresholds="best") # best thresholds
  best_thresholds = round(as.numeric(rownames(train_th$ci$specificity)),3)
  AUC = train_th$auc[1]
  train_thr_auc = data.frame(best_thresholds,AUC)
  
  train_thr_se_sp = data.frame(train_th$thresholds,
                               train_th$sensitivities,train_th$specificities) # sp、 se 、 threshold 
  colnames(train_thr_se_sp) = c('thresholds' ,'sensitivities ','specificities')
  merge_output_dat = list(train_thr_auc,train_thr_se_sp)
  txt_path = paste0(my_path,'/auc_thresholds.txt')
  capture.output(merge_output_dat, file = txt_path)
  # save img 
  png_auc = paste0(my_path,'/full-data-train_auc.png')
  png(png_auc,width = 1280,height = 720)
  rocobj <- plot.roc(train_y, pred_train_x,
                     
                     main="Confidence intervals of specificity/sensitivity", percent=TRUE,
                     
                     ci=TRUE, of="se", # ci of sensitivity
                     
                     specificities=seq(0, 100, 5), # on a select set of specificities
                     
                     ci.type="shape", ci.col="#1c61b6AA") # plot the CI as a blue shape
  
  plot(ci.sp(rocobj, sensitivities=seq(0, 100, 5)), # ci of specificity
       
       type="bars") # print this one as bars
  
  dev.off()
}

# predict new data-----------------------------------------------------------------------------------------------------
predict_function <-function(test_df,my_model,my_cutoff = 0.44){
  my_model = readRDS(my_model)# load model
  pred_score = as.data.frame(predict(my_model,test_df,type ='response'))
  pred_class = ifelse(pred_score$s0>my_cutoff,"Malignant","Benign")
  pred_score$pred_class = pred_class
  return(pred_score)
}

test_model <- function(data_input,label_input,best_alpha = 0.7,best_lambda=0.01){
  if(nchar(dirname(data_input))==1){
    my_path = getwd()
  }else{
    my_path = dirname(data_input)
  }
  #
  data_preoutput_output = data_pre(data_input,label_input) #  call the function : data_pre
  
  train_set = data_preoutput_output$train_set
  test_set = data_preoutput_output$test_set # return the values from data_pre
  
  train_x = as.matrix(subset(train_set,select = -label)) # train set
  train_y = train_set$label
  
  test_x = as.matrix(subset(test_set,select = -label))
  test_y = test_set$label
  
  fit_model = glmnet(train_x, train_y, family="binomial",alpha=0.3, lambda = 0.03397539) # model
  # save model
  my_model_path = paste0(my_path,"glmnet_model.rds")
  saveRDS(fit_model, file=my_model_path)
  pred_test_x= as.data.frame(predict(fit_model,test_x,type ='response'))$s0
  # save txt
  my_level = unique(train_y)
  test_th = roc(test_y, pred_test_x,levels = my_level,ci=TRUE, 
                of="thresholds",thresholds="best") # best thresholds
  best_thresholds = round(as.numeric(rownames(test_th$ci$specificity)),3)
  AUC = test_th$auc[1]
  train_thr_auc = data.frame(best_thresholds,AUC)
  
  train_thr_se_sp = data.frame(test_th$thresholds,
                               test_th$sensitivities,test_th$specificities) # sp、 se 、 threshold 
  colnames(train_thr_se_sp) = c('thresholds' ,'sensitivities ','specificities')
  merge_output_dat = list(train_thr_auc,train_thr_se_sp)
  txt_path = paste0(my_path,'/auc_thresholds.txt')
  capture.output(merge_output_dat, file = txt_path)
  # save img 
  png_auc = paste0(my_path,'/full-data-train_auc.png')
  png(png_auc,width = 1280,height = 720)
  rocobj <- plot.roc(test_y, pred_test_x,
                     
                     main="Confidence intervals of specificity/sensitivity", percent=TRUE,
                     
                     ci=TRUE, of="se", # ci of sensitivity
                     
                     specificities=seq(0, 100, 5), # on a select set of specificities
                     
                     ci.type="shape", ci.col="#1c61b6AA") # plot the CI as a blue shape
  
  plot(ci.sp(rocobj, sensitivities=seq(0, 100, 5)), # ci of specificity
       
       type="bars") # print this one as bars
  
  dev.off()
}


#input path
data_input_path = "/home/binyang_ni/R/R_data/glm_project_data/beta_binomial/LC002_tissue_145samples_betabinom_pv.tsv"
label_input_path = "/home/binyang_ni/R/R_data/glm_project_data/beta_binomial/label_beta_binomial.txt"
model_path = paste0(dirname(data_input_path),"/glmnet_model.rds")


select_par_output = select_par(data_input_path,label_input_path,n_times=5,my_step=0.1) # call select_par()
best_alpha = select_par_output$best_alpha
best_lambda = select_par_output$best_lambda

best_alpha = 0.4
best_lambda = 0.013
my_cutoff = 0.366

CV_threshold(data_input_path,label_input_path,best_alpha = best_alpha,best_lambda=best_lambda)# call CV_threshold()
train_model(data_input_path,label_input_path,best_alpha = best_alpha,best_lambda=best_lambda)# call train_model()
# get the test set
output_new = data_pre(data_input_path,label_input_path)
test_data = output_new$test_set
test_x = as.matrix(subset(test_data,select = -label)) # train set
test_y = test_data$label
# predict the test set
pred_score = predict_function(test_x,model_path,my_cutoff = 0.366)
pred_score$test_y = test_y
colnames(pred_score) = c("score","predict","TNM_stage")
table(test_y,pred_score$predict)
# ROC curve
rocobj <- plot.roc(test_y, pred_score$score,
                   
                   main="Confidence intervals of specificity/sensitivity", percent=TRUE,
                   
                   ci=TRUE, of="se", # ci of sensitivity
                   
                   specificities=seq(0, 100, 5), # on a select set of specificities
                   
                   ci.type="shape",print.auc=TRUE) # plot the CI as a blue shape

plot(ci.sp(rocobj, sensitivities=seq(0, 100, 5)), # ci of specificity
     
     type="bars") # print this one as bars


# my_model = readRDS("/home/binyang_ni/R/R_data/glm_project_data/glmnet_model.rds")# load model
# a11 = as.data.frame(as.matrix(my_model$beta))
# a11 = subset(a11,s0!=0)
# a22 = as.data.frame(my_model$a0)
# rownames(a22) = "a0"
# colnames(a22) = "s0"
# a33 = rbind(a22,a11)
# colnames(a33) = "beta"
# write.csv(a33,"/home/binyang_ni/R/R_data/glm_project_data/beta.csv")











