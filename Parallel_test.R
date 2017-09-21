# 并行计算euler14问题
# 自定义函数以返回原始数值和步数
func <- function(x) {
  n = 1
  raw <- x
  while (x > 1) {
    x <- ifelse(x%%2==0,x/2,3*x+1)
    n = n + 1
  }
  return(c(raw,n))
}

library(parallel)
# 用system.time来返回计算所需时间
timestart5<-Sys.time()
system.time({
  x <- 1:1e6
  cl <- makeCluster(4)  # 初始化四核心集群
  results <- parLapply(cl,x,func) # lapply的并行版本
  res.df <- do.call('rbind',results) # 整合结果
  stopCluster(cl) # 关闭集群
})
# 找到最大的步数对应的数字
res.df[which.max(res.df[,2]),1]
timeend6<-Sys.time()
print(timeend6-timestart5)



library(foreach)

timestart<-Sys.time()
####这块写你要运行的程序
# 非并行计算方式，类似于sapply函数的功能
x <- foreach(x=1:100000,.combine='rbind') %do% func(x)
timeend<-Sys.time()
print(timeend-timestart)
# 启用parallel作为foreach并行计算的后端

timestart3<-Sys.time()
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
# 并行计算方式
x <- foreach(x=1:100000,.combine='rbind') %dopar% func(x)
stopCluster(cl)
timeend4<-Sys.time()
print(timeend4-timestart3)



# 随机森林的并行计算
library(randomForest)
cl <- makeCluster(4)
registerDoParallel(cl)
rf <- foreach(ntree=rep(25000, 4), 
              .combine=combine,
              .packages='randomForest') %dopar%
  randomForest(Species~., data=iris, ntree=ntree)
stopCluster(cl)

rf2 = randomForest(Species~., data=iris, ntree=100000)
