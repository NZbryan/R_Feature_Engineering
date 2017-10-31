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
# index_9390 = input_dat$region
# input_dat_9390 = as.data.frame(input_dat)
input_dat_9390 = input_dat[-c(1),]
rownames(input_dat_9390) = input_dat_9390$Name
input_dat_9390 = input_dat_9390[,-c(1,2,3,4)]
### extract 1207 region from ' 9390 data '
# input_region1207 = fread('D:/R_workspace/LC006/LC002_tissue_beta_oddsratio_20170615.csv',check.names=FALSE,stringsAsFactors=FALSE)
model_path = '/home/binyang_ni/R/R_data/LC_model/LC002_tissue_Benign_IA_20170615.rds'
LC002_model = readRDS(model_path)# load model
a11 = as.data.frame(as.matrix(LC002_model$beta))
# a11 = subset(a11,s0!=0)

input_region1207 = rownames(a11)

### 1207 data unifNorm
data_1207 = input_dat_9390[input_region1207,]
LC006_index1207 = rownames(data_1207)
data_1207 = apply(data_1207,2,as.numeric)
data_1207 = 2^data_1207
data_1207 = data_1207-1

data_1207_unifNorm = apply(data_1207,2,unifNorm)
rownames(data_1207_unifNorm) = LC006_index1207

# write.csv(data_1207_unifNorm,
#           'D:/R_workspace/LC006/LC006_tissue_unifNorm_20170630.txt',sep='\t',quote = TRUE)

### test new data

test_x = t(data_1207_unifNorm)
test_x = test_x[,input_region1207]

pred_1= as.data.frame(predict(LC002_model,test_x,type ='response')) #   predict validation data
colnames(pred_1) = c("score")

pred_1$predict_class = ifelse(pred_1$score>  0.3,'Malignant','Benign')
# confusionMatrix(pred_1$predict_class,pred_1$TNMstage,positive ='Malignant')
write.csv(pred_1, paste0(my_path,'/score.csv'))

cat("output file: ", sprintf("%s", paste0(my_path,'/score.csv')), "\n sucess!\n")




