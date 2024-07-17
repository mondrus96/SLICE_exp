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

#### TWO THINGS TO EXPLORE ####
# 1: Normalization of outputs
# 2: Tune to similar sparsity and rank
# 3: Use CV once, set for all remaining

# Set params
params <- list(sli = c(0.25, 9),
               nnlvg = c(0.35, 0.014),
               rclvg = c(0.1, 9),
               tg = c(0.01, 0.2))

# Define a list of model functions
models <- list(
  sli = slice,
  nnlvg = nnlvg,
  rclvg = rclvg,
  tg = ebic.tg
)

# Run estimation of first iteration
set.seed(123)
X <- Xs[["Famous"]]
X <- X[sample(nrow(X)),] # Shuffle rows
mid <- nrow(X)/2
train <- 1:mid; test <- (mid + 1):nrow(X)

# Loop over models
model_outs <- list()
for(i in seq_along(models)){
  print(names(models)[i])
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
  print(paste0(names(models)[i], " - train s: ", 
               sum(train_s == 0)/length(train_s),
               ", test s: ", sum(test_s == 0)/length(test_s)))
  
  if(names(models)[i] != "tg"){
    train_r <- train_out$L; test_r <- test_out$L
    train_r <- sum(eigen(train_r)$values > 1e-3)
    test_r <- sum(eigen(test_r)$values > 1e-3)
    print(paste0(names(models)[i], " - train r: ", 
                 train_r, ", test r: ", test_r))
  }
}
    
# 100 iterations of 2-fold CV
for(cond in conds){ # Loop over conditions
  for(iter in 2:100){ # Loop over iterations
    set.seed(123*iter)
    X <- Xs[[cond]]
    X <- X[sample(nrow(X)),] # Shuffle rows
    mid <- nrow(X)/2
    train <- 1:mid; test <- (mid + 1):nrow(X)
    
    sli_train <- slice(cov(X[train,]), cvsli$rho, cvsli$r) # Train
    rclvg_train <- rclvg(cov(X[train,]), cvrclvg$rho, cvrclvg$r)
    nnlvg_train <- nnlvg(cov(X[train,]), cvnnlvg$rho, cvnnlvg$gamma)
    
    sli_test <- slice(cov(X[test,]), cvsli$rho, cvsli$r) # Test
    rclvg_test <- rclvg(cov(X[test,]), cvrclvg$rho, cvrclvg$r)
    nnlvg_test <- nnlvg(cov(X[test,]), cvnnlvg$rho, cvnnlvg$gamma)
    
    out1 <- frob.normalize(nnlvg_train$L)
    out2 <- frob.normalize(nnlvg_test$L)
    out3 <- frob.normalize(sli_train$L)
    out4 <- frob.normalize(sli_test$L)
    
    print(norm(out1 - out2, "F"))
    print(norm(out3 - out4, "F"))
    
    #sli1 <- eigen(sli_train$L)$vectors[,1]; sli2 <- eigen(sli_test$L)$vectors[,1]
    #nnlvg1 <- eigen(nnlvg_train$L)$vectors[,1]; nnlvg2 <- eigen(nnlvg_test$L)$vectors[,1]
    
    #sli_ang <- c(sli_ang, sintheta(sli1, sli2)); nnlvg_ang <- c(nnlvg_ang, sintheta(nnlvg1, nnlvg2))
    #sli_ang <- sintheta(sli1, sli2); nnlvg_ang <- sintheta(nnlvg1, nnlvg2)
    #print(paste0("slice:", mean(sli_ang), " nnvlg:", mean(nnlvg_ang)))
  }
}
