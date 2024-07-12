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
                         cond, "_hyp.txt"), ",", header = FALSE)
    X <- scale(X) # Scale data
    Xs[[cond]] <- rbind(Xs[[cond]], X)
  }
}

# Run estimation
set.seed(123)
cvslis <- slis <- list()
for(cond in conds){
  X <- Xs[[cond]]
  #cvslis[[cond]] <- cv.slice(X, rhos = logseq(0.05, 0.1, 4),
  #                           rs = 2:10, Sest = "huge_glasso")
  cvslis[[cond]] <- cv.slice(X, folds = 5, rhos = 0.2,
                             rs = 7:10, Sest = "huge_glasso")
  slis[[cond]] <- slice(cov(X), cvslis[[cond]]$rho, 
                        cvslis[[cond]]$r, Sest = "huge_glasso")
}
save(slis, cvslis, file = "ests.rda")