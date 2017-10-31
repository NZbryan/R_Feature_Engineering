#!/usr/bin/env Rscript
options(warn=-1) 
# pkgs = c('glmnet','caret','tidyr')
pkgs = c('glmnet')
choose_repos = 'http://mirrors.tuna.tsinghua.edu.cn/CRAN'
for(i in pkgs){
  if(!suppressWarnings(require(i,character.only = TRUE,warn.conflicts = F, quietly = T))){
    install.packages(i,repos=choose_repos)
    require(i,character.only = TRUE,warn.conflicts = F, quietly = T)
  }
}

args<-commandArgs(T)

if(nchar(dirname(args[1]))==1){
  my_path = getwd()
}else{
  my_path = dirname(args[1]) # get the file path
}

### unifNorm function
unifNorm<- function(x, trunPoint=0){
  require(Hmisc)
  idx=x>trunPoint
  res = ecdf(x[idx])
  y=x
  y[idx]=res(x[idx])
  return(y)
}
set.seed(20170415)
### 9390 data
input_dat = read.csv(args[1],
                     sep = '\t',check.names=FALSE,stringsAsFactors=FALSE)

input_dat_9390 = input_dat[-c(1),]
index_9390 = input_dat_9390$Name
rownames(input_dat_9390) = index_9390
input_dat_9390 = input_dat_9390[,-c(1,2,3,4)]

### 9390 data unifNorm
data_9390_CPM= apply(input_dat_9390,2,as.numeric)
data_9390_CPM = 2^data_9390_CPM
data_9390_CPM = data_9390_CPM-1

data_9390_unifNorm = apply(data_9390_CPM,2,unifNorm)
rownames(data_9390_unifNorm) = index_9390

### load 
model_path = '/home/binyang_ni/R/R_data/LC_model/LC005_train_allsample_20170703.rds'
LC002_model = readRDS(model_path)# load model
a11 = as.data.frame(as.matrix(LC002_model$beta))

### test new data
test_x = t(data_9390_unifNorm)
test_x = test_x[,rownames(a11)]

pred_1= as.data.frame(predict(LC002_model,test_x,type ='response')) #   predict validation data
colnames(pred_1) = c("score")

pred_1$predict_class = ifelse(pred_1$score>  0.525,'Malignant','Benign')
# confusionMatrix(pred_1$predict_class,pred_1$TNMstage,positive ='Malignant')
write.csv(pred_1, paste0(my_path,'/LC005predict_score.csv'))

cat("output file: ", sprintf("%s", paste0(my_path,'/LC005predict_score.csv')), " sucess!\n")

