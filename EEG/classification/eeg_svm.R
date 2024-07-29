library(e1071)
library(MASS)
library(parallel)
sapply((paste0("../../Core/", list.files("../../Core/"))), source)

### Group-wise ###
# Load data
Xs <- list()
conds <- c("Famous", "Scrambled", "Unfamiliar")
dirs <- list.dirs("..", recursive = FALSE)
dirs <- dirs[!grepl("../classification", dirs)]
# Loop over conditions
for(cond in conds){
  # Loop over subjects
  for(dir in dirs){
    X <- read.csv(paste0(dir, "/meeg/", 
                         cond, ".txt"), ",", header = FALSE)
    Xs[[cond]] <- rbind(Xs[[cond]], X)
  }
  Xs[[cond]] <- scale(Xs[[cond]]) # Scale data
}

# Define a function for the inner loop
inner.loop <- function(iter, model, train_cov, test_cov, n2, param1, param2){
  # Train and test covariances
  X_train_ <- cov(mvrnorm(n2, rep(0, ncol(train_cov)), train_cov))
  X_test_ <- cov(mvrnorm(n2, rep(0, ncol(test_cov)), test_cov))
  
  # Define a named list of functions for each model
  model_funcs <- list(
    slice = slice,
    nnlvg = nnlvg,
    rclvg = rclvg
  )
  
  # Train and test data
  train <- model_funcs[[model]](X_train_, param1, param2)
  test <- model_funcs[[model]](X_test_, param1, param2)
  
  # Features
  X_train_row <- eigen(train$L)$vectors[,1]
  X_test_row <- eigen(test$L)$vectors[,1]
  
  list(X_train_row = X_train_row, X_test_row = X_test_row)
}

# Parameter generation
n2 <- 200
num_cores <- detectCores() # Use all cores
cl <- makeCluster(num_cores)
clusterExport(cl, c("mvrnorm", "slice", "eigs", "objective",
                    "inner.loop", "glasso", "n2", "cov", 
                    "eigen", "svds", "nucl_norm", "logL", 
                    "makePD", "isPD", "nucl_shr", "Diagonal", 
                    "nnlvg", "L1_shr", "rclvg", "Estep",
                    "Mstep", "forceSymmetric"))
clusterEvalQ(cl, {
  library(MASS)
  library(e1071)
  library(Matrix)
})

reps <- 100
lambdas <- logseq(1e-4, 1e-8, 5)
rs <- 6:10
param_len <- length(lambdas)

mats <- list()
models <- c("slice", "nnlvg", "rclvg")

# Loop over models
for(model in models){
  mat <- matrix(NA, param_len, param_len)
  for(i in 1:param_len){
    for(j in 1:param_len){
      print(paste0("param1 iter: ", i, ", param2 iter: ", j))
      
      result <- tryCatch({
        set.seed(123) # For sampling
        X_train <- y_train <- X_test <- y_test <- c()
        for(cond in conds){
          X <- Xs[[cond]]
          X <- X[sample(nrow(X)), ]
          n1 <- nrow(X)/2
          train_cov <- cov(X[1:n1, ]) # Covariances used
          test_cov <- cov(X[(n1 + 1):nrow(X), ])
          
          if(model == "nnlvg"){
            param1 <- lambdas[i]; param2 <- lambdas[j]
          } else{
            param1 <- lambdas[i]; param2 <- rs[j]
          }
          
          clusterSetRNGStream(cl, 123) # For parallel processes
          results <- parLapply(cl, 1:reps, inner.loop, model,
                               train_cov, test_cov, n2, param1, param2)
          
          # Collect results
          X_train_rows <- lapply(results, function(x) x$X_train_row)
          X_test_rows <- lapply(results, function(x) x$X_test_row)
          
          X_train <- rbind(X_train, do.call(rbind, X_train_rows))
          X_test <- rbind(X_test, do.call(rbind, X_test_rows))
          
          y_train <- c(y_train, rep(cond, reps))
          y_test <- c(y_test, rep(cond, reps))
        }
        
        # Shuffle data
        row_perm <- sample(nrow(X_train))
        X_train <- X_train[row_perm,]; y_train <- y_train[row_perm]
        row_perm <- sample(nrow(X_test))
        X_test <- X_test[row_perm,]; y_test <- y_test[row_perm]
        
        # SVM classifier
        svm_model <- svm(X_train, as.factor(y_train), 
                     type = "C-classification", kernel = "radial")
        y_pred <- as.character(predict(svm_model, X_test))
        accuracy <- mean(y_pred == y_test)
        print(accuracy)
      }, error = function(e) {
        print(paste("Error occurred:", conditionMessage(e)))
        accuracy <- NA
      })
      
      # Store the result in the matrix
      mat[i, j] <- result
    }
  }
  mats[[model]] <- mat
  save(mats, file = "svmmats.rda")
}
stopCluster(cl)