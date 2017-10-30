library(glmnet)
library(caret)
library(data.table)
options(warn=-1) 
#
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
  
  lable_info_train = subset(lable_info, lable_info$IsTraining==1)
  lable_info_test = subset(lable_info, lable_info$IsTraining==0)
  
  train_set = input_dat[lable_info_train$SampleID,]
  test_set = input_dat[lable_info_test$SampleID,]
  
  train_set$label = lable_info_train$label
  test_set$label = lable_info_test$label
  return(list(train_set=train_set,test_set=test_set ))
}

# CV
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
      fit_bi = glmnet(train_x, train_y, family="binomial",alpha=best_alpha, lambda = best_lambda) # model
      pred_vali= as.data.frame(predict(fit_bi,vali_x,type ='response'))$s0 #   predict validation data
      
      my_level = unique(train_y)
      vali_th = roc(vali_y, pred_vali,levels = my_level,ci=TRUE, of="thresholds",thresholds="best") # best thresholds
      best_thresholds = round(as.numeric(rownames(vali_th$ci$specificity)),3)
      AUC = vali_th$auc[1]
      vali_thr_auc = data.frame(best_thresholds,AUC)
      my_underline = paste0('-------------------------------',k-1,i,'-----------------------------------')
      all_thresholds = append(all_thresholds,best_thresholds)
      
      vali_thr_se_sp = data.frame(vali_th$thresholds,vali_th$sensitivities,vali_th$specificities) # sp、 se 、 threshold 
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


# CV_threshold('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt','/home/binyang_ni/R/R_data/glm_project_data/labe_info.txt')


