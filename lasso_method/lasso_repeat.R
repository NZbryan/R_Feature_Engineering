library(glmnet)

feature_selection <- function(train_x,train_y,sampling_times = 50){
  beta_cal = data.frame(feature_name=colnames(train_x))
  for(i in seq(sampling_times)){
    set_myseed = as.integer(paste(201701017,i,sep = ""))
    set.seed(set_myseed)
    # sampling
    sub1 =  sample(nrow(train_x), as.integer(nrow(train_x)*0.75), replace = F)
    train_x = train_x[sub1,]
    train_y = train_y[sub1]
    # cv
    fit_cv = cv.glmnet(train_x, train_y, family="binomial",alpha=1,nfolds=10) # run the model
    # fit model
    fit_model = glmnet(train_x, train_y, family="binomial",
                     alpha=1, lambda = fit_cv$lambda.1se) # model
    
    feature_important = as.data.frame(as.matrix(fit_model$beta))
    feature_important = feature_important[colnames(train_x),]
    beta_cal = cbind(beta_cal,feature_important)
    colnames(beta_cal)[i+1] = paste('time_',i,sep = "")

  }
  return(beta_cal)
}


input_path = '/home/binyang_ni/R/R_data/jiacheng_chuan/201709/input_dat_scale_201709.csv'
label_path = '/home/binyang_ni/R/R_data/jiacheng_chuan/201709/label_full.tsv'
data_pre <- function(data_input,label_input){
  # input_dat = fread('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt',stringsAsFactors=FALSE, check.names=TRUE ,sep = '\t')
  # input_dat = read.csv('/home/binyang_ni/python/data/LC002_Tissue_QCpassed_093Up_20170309.CPM.txt',header = T,
  #                      row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t')
  # lable_info = fread('/home/binyang_ni/R/R_data/glm_project_data/labe_info.txt',stringsAsFactors=FALSE, check.names=TRUE)
  input_dat = read.csv(data_input,header = T,
                       row.names = 1,stringsAsFactors = FALSE,check.names = FALSE,sep = '\t')
  # input_dat = read.csv(data_input,header = T,
  #                      row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
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
data_out = data_pre(input_path,label_path)
dat = data_out$train_set
train_x = as.matrix(subset(dat,select = -label))
train_y = dat$label

output = feature_selection(train_x,train_y,sampling_times = 500)

