#!/usr/bin/env Rscript

# rm(list = ls(all = TRUE))
options(warn=-1)
options(digits = 4)
### load the library
pkgs = c('ROCR','argparse')
choose_repos = 'http://mirrors.tuna.tsinghua.edu.cn/CRAN'
for(i in pkgs){
  if(!suppressWarnings(require(i,character.only = TRUE,warn.conflicts = F, quietly = T))){
    install.packages(i,repos=choose_repos)
    require(i,character.only = TRUE,warn.conflicts = F, quietly = T)
  }
}

parser <- ArgumentParser(description='Feature Selection---forward greedy algorithm')

parser$add_argument('-train', dest='train_data_arg', metavar='FILE',required=TRUE,
                    help='train data file, the tab-deliminated data matrix file,row is sample and  column is feature ')

parser$add_argument('-test', dest='test_data_arg', metavar='FILE',required=TRUE, 
                    help='test data file, the tab-deliminated data matrix file,row is sample and  column is feature ')

parser$add_argument('-train_threshold',dest='train_thre',action='store',type = 'double',
                    default = 0.95,help='The AUC threshold of train data')

parser$add_argument('-test_threshold',dest='test_thre',action='store',type = 'double',
                    default = 0.85,help='The AUC threshold of test data')

parser$add_argument('-o', dest='outpath', metavar='path', type='character',
                    help='output path')

args <- parser$parse_args()

train_data = read.csv(args$train_data_arg,sep = '\t')
train_y = factor(train_data$label)
train_x = subset(train_data,select = -label)

test_data = read.csv(args$test_data_arg,sep = '\t')
test_y = factor(test_data$label)
test_x = subset(test_data,select = -label)


Calc_auc = function(score,label){
  label = as.character(label)
  pred <- prediction(score, label)
  AUC = performance(pred, "auc")@y.values[[1]]
  if(AUC<0.5){
    new_label = ifelse(label==unique(label)[1],unique(label)[2],unique(label)[1])
    pred <- prediction(score, new_label)
    AUC = performance(pred, "auc")@y.values[[1]]
  }
  return(AUC)
}

getmax = function(feature_name,AUC_list){
  AUC_max=max(AUC_list,na.rm=TRUE)
  idmax=AUC_list>=AUC_max
  feature_selected=feature_name[idmax][1]
  return(list(feature_selected=feature_selected,AUC_max=AUC_max))
}


sga = as.numeric(lapply(seq(ncol(train_x)),function(i) Calc_auc(train_x[,i],train_y)))
single_feature_auc = getmax(colnames(train_x),sga)
first_f = single_feature_auc$feature_selected
first_auc = single_feature_auc$AUC_max
all_feature = setdiff(colnames(train_x),first_f)



basic_f = first_f
basic_auc = first_auc
j=0
test_auc_first = Calc_auc(test_x[,first_f],test_y)
auc_list_test = test_auc_first
MAX_AUC_FUN_train_test = function(loop_feature){
  j<<-j+1
  auc_list_train = c()
  for(i in seq(length(loop_feature))){
    model_f = append(basic_f,loop_feature[i])
    # fit_cv = cv.glmnet(train_x[,model_f], train_y, family="binomial",alpha=1,nfolds=10,type.measure="auc") # run the model
    # train_predict = predict(fit_cv, train_x[,model_f], s=fit_cv$lambda.min)
    train_daf = train_x[,model_f]
    train_daf$label = train_y
    fit_glm = glm(label ~.,family=binomial(link='logit'),data=train_daf)
    train_predict =  predict(fit_glm,newdata=train_x[,model_f],type='response')
    
    train_auc = Calc_auc(train_predict,train_y)
    auc_list_train = append(auc_list_train,train_auc)
    
  }
  
  get_max = getmax(loop_feature,auc_list_train)
  AucMax = get_max$AUC_max
  AucMax_f = get_max$feature_selected
  
  basic_f <<- append(basic_f,AucMax_f)
  basic_auc <<- append(basic_auc,AucMax)
  loop_feature = setdiff(loop_feature,AucMax_f)
  
  # test
  train_daf = train_x[,basic_f]
  train_daf$label = train_y
  fit_glm = glm(label ~.,family=binomial(link='logit'),data=train_daf)
  
  test_daf = test_x[,basic_f]
  test_daf$label = test_y
  test_predict = predict(fit_glm,newdata=test_x[,basic_f],type='response')
  test_auc = Calc_auc(test_predict,test_y)
  auc_list_test <<- append(auc_list_test,test_auc)
  
  cat(sprintf("feature个数: %d, train auc: %f, test auc: %s\n",j+1,AucMax,test_auc))
  
  # outp_df <<- cbind()
  if(AucMax>args$train_thre & test_auc>args$test_thre){
    return(data.frame(feature_selected=basic_f,AccumTrain_AUC=basic_auc,AccumTest_AUC=auc_list_test))
  }else{
    # return(MAX_AUC_FUN_train(loop_feature))
    MAX_AUC_FUN_train_test(loop_feature)
  }
}

cat('-----------------------output----------------------------------------------\n')

first_test_auc = Calc_auc(test_x[,first_f],test_y)
cat(sprintf("feature个数: %d, train auc: %f, test auc: %s\n",j+1,first_auc,first_test_auc))

outp = MAX_AUC_FUN_train_test(all_feature)

if(is.null(args$outpath)){
  my_path = getwd()
}else{
  dir.exists(args$outpath) || dir.create(args$outpath)
  my_path = args$outpath
}

write.csv(outp,paste0(my_path,'/output.csv'))


