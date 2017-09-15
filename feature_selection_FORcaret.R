library(caret)
data(BloodBrain)
set.seed(1)
RFwithGAM <- sbf(bbbDescr, logBBB,
                 sbfControl = sbfControl(functions = rfSBF,
                                         verbose = FALSE,
                                         seeds = sample.int(100000, 11),
                                         method = "cv"))
RFwithGAM



if(!suppressWarnings(require(C50))){
  install.packages('C50')
  require(C50)
}


sbfControls_nb <- sbfControl(
  functions = nbSBF,
  method = 'boot')

fs_nb <- sbf(x = churnTrain[,-20],
             y = churnTrain[,20],
             sbfControl = sbfControls_nb)



###
sbfControls_rf <- sbfControl(
  functions = rfSBF,
  method = 'cv',
  repeats = 5)

fs_rf <- sbf(x = churnTrain[,-20],
             y = churnTrain[,20],
             sbfControl = sbfControls_rf)


rfeControls_nb <- rfeControl(
  functions = nbFuncs,
  method = 'boot')

fs_nb2 <- rfe(x = churnTrain[,-20],
             y = churnTrain[,20],
             sizes = seq(4,19,2),
             rfeControl = rfeControls_nb)



#构建rfe函数的控制参数(使用随机森林函数和10重交叉验证抽样方法，并抽取5组样本)
rfeControls_rf <- rfeControl(
  
  functions = rfFuncs,
  
  method = 'cv',
  
  repeats = 5)


#使用rfe函数进行特征选择

fs_nb3 <- rfe(x = churnTrain[,-20],
             
             y = churnTrain[,20],
             
             sizes = seq(4,19,2),
             
             rfeControl = rfeControls_rf)






library(caret)
library(mlbench)
library(Hmisc)
library(randomForest)

n <- 100
p <- 40
sigma <- 1
set.seed(1)

sim <- mlbench.friedman1(n, sd = sigma)
colnames(sim$x) <- c(paste("real", 1:5, sep = ""),
                     paste("bogus", 1:5, sep = ""))
bogus <- matrix(rnorm(n * p), nrow = n)
colnames(bogus) <- paste("bogus", 5+(1:ncol(bogus)), sep = "")
x <- cbind(sim$x, bogus)
y <- sim$y



normalization <- preProcess(x)
x <- predict(normalization, x)
x <- as.data.frame(x)
subsets <- c(1:5, 10, 15, 20, 25)



set.seed(10)

ctrl <- rfeControl(functions = lmFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(x, y,
                 sizes = subsets,
                 rfeControl = ctrl)

lmProfile


trellis.par.set(caretTheme())
plot(lmProfile, type = c("g", "o"))






rfRFE <-  list(summary = defaultSummary,
               fit = function(x, y, first, last, ...){
                 library(randomForest)
                 randomForest(x, y, importance = first, ...)
               },
               pred = function(object, x)  predict(object, x),
               rank = function(object, x, y) {
                 vimp <- varImp(object)
                 vimp <- vimp[order(vimp$Overall,decreasing = TRUE),,drop = FALSE]
                 vimp$var <- rownames(vimp)                  
                 vimp
               },
               selectSize = pickSizeBest,
               selectVar = pickVars)



if(suppressWarnings(require(C5sa0))){
  print('test')
}
