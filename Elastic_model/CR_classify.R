#!/usr/bin/env Rscript

pkgs = c('glmnet','caret','data.table','pROC','argparse')
choose_repos = 'http://mirrors.tuna.tsinghua.edu.cn/CRAN'
for(i in pkgs){
  if(!suppressWarnings(require(i,character.only = TRUE,warn.conflicts = F, quietly = T))){
    install.packages(i,repos=choose_repos)
    require(i,character.only = TRUE,warn.conflicts = F, quietly = T)
  }
}

parser <- ArgumentParser(description='CR003 tissue classification')

parser$add_argument('-s', '--sample', dest='sample', metavar='FILE',required=TRUE, nargs='+',
                    help='data file, the tab-deliminated data matrix file,row is region and  column is sample ')

parser$add_argument('-l', '--label', dest='label', metavar='FILE',required=TRUE, nargs='+',
                    help='label file, the tab-deliminated patient info. data')

parser$add_argument('-o', '--outdir', dest='outpath', metavar='path', type='character',
                    help=sprintf('output path [StepWiseRemoval.%s]', format(Sys.Date(), "%y%m%d")))

args <- parser$parse_args()


data_pre <- function(data_input,label_input){

  # raw data
  input_dat = read.csv(data_input[1],header = T,
                       row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t')
  if (length(data_input) > 1) {
    for (i in 2:length(data_input)) {
      input_dat = cbind(input_dat, read.csv(data_input[1],header = T,
                                     row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t'))
    }
  }
  
  input_dat = as.data.frame(t(input_dat))
  normalization = preProcess(input_dat)
  input_dat = predict(normalization, input_dat)
  input_dat = as.data.frame(input_dat)
  
  
  # label
  lable_info = fread(label_input[1],stringsAsFactors=FALSE, check.names=TRUE)
  if (length(label_input) > 1) {
    for (i in 2:length(label_input)) {
      lable_info = rbind(lable_info, fread(label_input[1],
                                          stringsAsFactors=FALSE, check.names=TRUE))
    }
  }
  
  
  colnames(lable_info) = c('SampleID','IsTraining','label')

  lable_info_train = subset(lable_info, lable_info$IsTraining==1)
  lable_info_test = subset(lable_info, lable_info$IsTraining==0)
  
  train_set = input_dat[lable_info_train$SampleID,]
  test_set = input_dat[lable_info_test$SampleID,]
  
  train_set$label = lable_info_train$label
  test_set$label = lable_info_test$label
  return(list(train_set=train_set,test_set=test_set,mean=normalization$mean,sd=normalization$std))
}

### parameter tuning
select_par <- function(data_input,label_input, output_path,n_times = 50,my_step = 0.01){
  
  # if(nchar(dirname(data_input))==1){
  #   my_path = getwd()
  # }else{
  #   my_path = dirname(data_input) # get the file path
  # }
  my_path = output_path
  data_preoutput_output = data_pre(data_input,label_input) #  call the function : data_pre
  df = data_preoutput_output$train_set # return the values from data_pre
  set.seed(201703233)
  
  train_x = as.matrix(subset(df,select = -label)) # train set
  train_y = df$label
  
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
              output=output)) #  return alpha閵?? lambda 閵?? and other relative values
}

### cross validation based on training data
CV_threshold <- function(data_input,label_input,output_path,best_alpha = 0.7,best_lambda=0.01){
  # if(nchar(dirname(data_input))==1){
  #   my_path = getwd()
  # }else{
  #   my_path = dirname(data_input)
  # }
  #
  my_path = output_path
  
  data_preoutput_output = data_pre(data_input,label_input) #data input
  train_set = data_preoutput_output$train_set
  
  pdf_path = paste0(my_path,'/auc_plot.pdf') # pdf path
  pdf(pdf_path,width =8.5 , height = 20)
  all_thresholds = c()
  auc_threshold_list =list()
  my_step = 9
  loop_count = 0
  sp = c()
  se = c()
  best_cutoff = c()
  auc_model = c()
  for(i in 0:my_step){
    set_myseed = as.integer(paste(20170323,i,sep = ""))
    set.seed(set_myseed)
    cv_Folds = createFolds(1:dim(train_set)[1],k=3)
    par(mfrow = c(3,1))
    
    for(k in 1:length(cv_Folds))
    { 
      loop_count = loop_count+1
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
      # save model
      my_model_path = paste0(my_path,"/glmnet_model",loop_count,".rds")
      saveRDS(fit_bi, file=my_model_path)
      my_level = unique(train_y)
      vali_th = roc(vali_y, pred_vali,levels = my_level,ci=TRUE, 
                    of="thresholds",thresholds="best") # best thresholds
      best_thresholds = round(as.numeric(rownames(vali_th$ci$specificity)),3)
      AUC = vali_th$auc[1]
      se = append(se,vali_th$ci$sensitivity[2])
      sp = append(sp,vali_th$ci$specificity[2])
      best_cutoff = append(best_cutoff,best_thresholds)
      auc_model = append(auc_model,AUC)
      vali_thr_auc = data.frame(best_thresholds,AUC)
      my_underline = paste0('-------------------------------',k-1,i,'-----------------------------------')
      all_thresholds = append(all_thresholds,best_thresholds)
      
      vali_thr_se_sp = data.frame(vali_th$thresholds,
                                  vali_th$sensitivities,vali_th$specificities) # sp銆? se 銆? threshold 
      colnames(vali_thr_se_sp) = c('thresholds' ,'sensitivities ','specificities')
      mergeouy_output_dat = list(my_underline,vali_thr_auc,vali_thr_se_sp)
      auc_threshold_list = append(auc_threshold_list,mergeouy_output_dat) # list the output
      # txt_path = paste0(my_path,'/auc_thresholds.txt')
      # capture.output(out_se_sp_cutoff_auc, file = txt_path,append = TRUE)
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
  out_se_sp_cutoff_auc = data.frame(se,sp,best_cutoff,auc_model)
  txt_path = paste0(my_path,'/out_se_sp_cutoff_auc.txt')
  capture.output(out_se_sp_cutoff_auc, file = txt_path)
}

### boostrap 3 times test model
boostrap_testModel <- function(data_input,label_input,output_path,best_alpha = 0.7,best_lambda=0.01){
  # output path
  my_path = output_path
  #
  data_preoutput_output = data_pre(data_input,label_input) #  call the function : data_pre
  train_set = data_preoutput_output$train_set # return the values from data_pre
  
  train_x = as.matrix(subset(train_set,select = -label)) # train set
  train_y = train_set$label
  # random sampling function
  random_sampling <- function(train_x,train_y,best_alpha,best_lambda,n){
    random_sampling_data = data.frame(sampleID=rownames(train_x),TNMstage=train_y)
    for(i in seq(n)){
      random_sampling_data[paste0('random_sampling_',i)] = NA
      set_myseed = as.integer(paste(20170921,i,sep = ""))
      set.seed(set_myseed)
      sub1 =  sample(nrow(train_x), as.integer(nrow(train_x)/3), replace = F)
      random_sampling_data[paste0('random_sampling_',i)][sub1,] = 0
      random_sampling_data[paste0('random_sampling_',i)][-sub1,] = 1
      # train 
      X_train = train_x[sub1,]
      y_train = train_y[sub1,]
      
      X_test = train_x[-sub1,]
      y_test = train_y[-sub1,]
      fit_model = glmnet(X_train, y_train, family="binomial",alpha=best_alpha, lambda = best_lambda) # model
      pred= as.data.frame(predict(fit_model,X_test,type ='response')) #   predict validation data
      pred$TNMstage = y_test
      colnames(pred) = c("score","TNMstage")
      # save score
      write.csv(pred, sprintf("%s/boostrap_score.xlsx", my_path), 
                row.names=F, sep="\t", quote=F,sheetName=paste0('score_random_sampling_',i),append = T)
      # save model
      saveRDS(fit_model,paste0(my_path,'/random_sampling_',i,'model.rds'))
      
      # ROC curve
      png_path = paste0(my_path,'/ROC_',i,'.png')
      png(png_path,width = 436,height = 438)
      rocobj <- plot.roc(y_test, pred$score,
                         
                         main="Confidence intervals of specificity/sensitivity", percent=TRUE,
                         
                         ci=TRUE, of="se", # ci of sensitivity
                         
                         specificities=seq(0, 100, 5), # on a select set of specificities
                         
                         ci.type="shape", ci.col="#1c61b6AA",print.auc=TRUE) # plot the CI as a blue shape
      
      plot(ci.sp(rocobj, sensitivities=seq(0, 100, 5)), # ci of specificity
           
           type="bars") # print this one as bars
      
      dev.off()
      
    }
    write.table(random_sampling_data, sprintf("%s/boostrap_score.xlsx", my_path), 
                row.names=F, sep="\t", quote=F,sheetName="random_sampling_total_info",append = T)
    
  }
  
  # random_sampling_data
  # random_sampling_data = random_sampling(train_x,train_y,3)
  # write.table(random_sampling_data, sprintf("%s/random_sampling_data.xlsx", my_path), row.names=F, sep="\t", quote=F,sheetName="random_sampling")
  random_sampling(train_x,train_y,3)
}


train_model <- function(data_input,label_input,output_path,best_alpha = 0.7,best_lambda=0.01){
  # output path
  my_path = output_path
  #
  data_preoutput_output = data_pre(data_input,label_input) #  call the function : data_pre
  train_set = data_preoutput_output$train_set # return the values from data_pre
  
  train_x = as.matrix(subset(train_set,select = -label)) # train set
  train_y = train_set$label
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




outp = select_par(args$sample,args$label,args$outpath,n_times = 30,my_step = 0.1)

CV_threshold(args$sample,args$label,args$outpath,best_alpha = outp$best_alpha,best_lambda=outp$best_lambda)






