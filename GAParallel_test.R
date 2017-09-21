
library(mlbench)
library(caret)
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


### GA
registerDoParallel(4) # Registrer a parallel backend for train
getDoParWorkers() # check that there are 4 workers

ga_ctrl <- gafsControl(functions = rfGA,
                       method = "repeatedcv",
                       repeats = 5,
                       genParallel=TRUE, # Use parallel programming
                       allowParallel = TRUE)

## Use the same random number seed as the RFE process
## so that the same CV folds are used for the external
## resampling. 
set.seed(10)
system.time(rf_ga <- gafs(x = x, y = y,
                          iters = 200,
                          gafsControl = ga_ctrl))




registerDoParallel(4) # Registrer a parallel backend for train
getDoParWorkers() # check that there are 4 workers

ga_ctrl <- gafsControl(functions = rfGA, # Assess fitness with RF
                       method = "cv",    # 10 fold cross validation
                       genParallel=TRUE, # Use parallel programming
                       allowParallel = TRUE)
## 
set.seed(10)
lev <- c("PS","WS")     # Set the levels

system.time(rf_ga3 <- gafs(x = trainX, y = y,
                           iters = 100, # 100 generations of algorithm
                           popSize = 20, # population size for each generation
                           levels = lev,
                           gafsControl = ga_ctrl))

