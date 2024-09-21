library(dplyr)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Load data
Xs <- list()
conds <- c("Famous", "Scrambled", "Unfamiliar")
dirs <- list.dirs(recursive = FALSE)
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

# Function for min max normalization
minmax.normalize <- function(mat) {
  min_val <- min(mat)
  max_val <- max(mat)
  mat <- (mat - min_val) / (max_val - min_val)
  return(mat)
}

frob.normalize <- function(mat) {
  return(mat/norm(mat, "F"))
}

# Set params
params <- list(sli = c(0.15, 9),
               nnlvg = c(0.1, 0.007),
               rclvg = c(0.15, 9),
               tg = c(0.15, 0.2))

# Define a list of model functions
models <- list(
  sli = slice,
  nnlvg = nnlvg,
  rclvg = rclvg,
  tg = ebic.tg
)

# Loop over models
X <- Xs[["Famous"]]
all_s <- all_r <- model_name <- ARI_score <- F1_score <- df <- c()
for(iter in 1:100){ # Loop over iterations
  print(paste0("ITER: ", iter))
  set.seed(123*iter)
  X <- X[sample(nrow(X)),] # Shuffle rows
  mid <- nrow(X)/2
  train <- 1:mid; test <- (mid + 1):nrow(X)
  #for(i in seq_along(models)){
  for(i in 1:2){
    model_name <- c(model_name, names(models)[i])
    model <- models[[i]]
    
    # Train and test
    train_out <- model(cov(X[train,]), 
                       params[[i]][1], 
                       params[[i]][2])
    test_out <- model(cov(X[test,]), 
                      params[[i]][1], 
                      params[[i]][2])
    model_outs[[names(models)[i]]] <- list(train = train_out, 
                                           test = test_out)
    
    # Compare sparsity of S and rank of L
    train_s <- train_out$S; test_s <- test_out$S
    train_s <- train_s[upper.tri(train_s)]; test_s <- test_s[upper.tri(test_s)]
    avg_s <- (sum(train_s == 0)/length(train_s) + 
                sum(test_s == 0)/length(test_s))/2
    all_s <- c(all_s, avg_s)
    
    if(names(models)[i] != "tg"){
      train_r <- train_out$L; test_r <- test_out$L
      train_r <- sum(eigen(train_r)$values > 1e-3)
      test_r <- sum(eigen(test_r)$values > 1e-3)
      avg_r <- (train_r + test_r)/2
      all_r <- c(all_r, avg_r)
    } else{
      all_r <- c(all_r, NA)
    }
    
    # Compare train and test
    F1_score <- c(F1_score, 
                  F1score(train_out$S, test_out$S))
    if(names(models[i]) == "tg"){
      ARI_score <- c(ARI_score, NA)
    } else{
      train_eigL <- eigen(train_out$L)$vectors[,1:train_r]
      test_eigL <- eigen(test_out$L)$vectors[,1:test_r]
      train_c <- kmeans(train_eigL, train_r, 10000)$cluster
      test_c <- kmeans(test_eigL, test_r, 10000)$cluster
      ARI_score <- c(ARI_score, ari(train_c, test_c))
    }
    df <- rbind(df, data.frame(model_name, all_s, all_r, F1_score, ARI_score))
  }
  mean_df <- df %>%
    group_by(model_name) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  print(mean_df)
}
