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
  
  # Train and test data
  train <- slice(X_train_, param1, param2)
  test <- slice(X_test_, param1, param2)
  
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
                    "nnlvg", "L1_shr"))
clusterEvalQ(cl, {
  library(MASS)
  library(e1071)
  library(Matrix)
})

reps <- 100
lambdas <- logseq(0.1, 1e-5, 5)
rs <- 6:10
model <- "slice"
mat <- matrix(NA, length(lambdas), length(rs))
colnames(mat) <- format(lambdas, scientific = TRUE)
rownames(mat) <- format(rs, scientific = TRUE)

# Loop over models
for(i in 1:length(rs)){
  for(j in 1:length(lambdas)){
    set.seed(123) # For sampling
    X_train <- y_train <- X_test <- y_test <- c()
    for(cond in conds){
      X <- Xs[[cond]]
      X <- X[sample(nrow(X)), ]
      n1 <- nrow(X)/2
      train_cov <- cov(X[1:n1, ]) # Covariances used
      test_cov <- cov(X[(n1 + 1):nrow(X), ])
      
      clusterSetRNGStream(cl, 123) # For parallel processes
      results <- parLapply(cl, 1:reps, inner.loop, model,
                           train_cov, test_cov, n2, 
                           lambdas[j], rs[i])
      
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
    accuracy <- mean(y_pred == y_test)*100
    print(accuracy)
    
    # Store the result in the matrix
    mat[i, j] <- accuracy
  }
}
save(mat, file = "svmmat.rda")

# For getting best results from all combinations
best <- which(mat == max(mat), arr.ind = TRUE)

set.seed(123) # For sampling
X_train <- y_train <- X_test <- y_test <- c()
for(cond in conds){
  X <- Xs[[cond]]
  X <- X[sample(nrow(X)), ]
  n1 <- nrow(X)/2
  train_cov <- cov(X[1:n1, ]) # Covariances used
  test_cov <- cov(X[(n1 + 1):nrow(X), ])
  
  clusterSetRNGStream(cl, 123) # For parallel processes
  results <- parLapply(cl, 1:reps, inner.loop, model,
                       train_cov, test_cov, n2, lambdas[best[2]], rs[best[1]])
  
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

# Append to svmpreds
svmpreds <- list(y_pred = y_pred, y_test = y_test)
save(svmpreds, file = "svmpreds.rda")
  
# STOP CLUSTER
stopCluster(cl)