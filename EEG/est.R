sapply((paste0("../Core/", list.files("../Core/"))), source)

### Group-wise ###
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

# Run estimation
set.seed(123)
cvslis <- slis <- cvnnlvgs <- nnlvgs <- cvrclvgs <- rclvgs <- tgs <- list()
for(cond in conds){
  X <- Xs[[cond]]
  
  # SLICE
  cvslis[[cond]] <- cv.slice(X, rhos = logseq(0.01, 0.5, 5), rs = 2:10)
  slis[[cond]] <- slice(cov(X), cvslis[[cond]]$rho,
                        cvslis[[cond]]$r)
  
  # nnLVGLASSO
  cvnnlvgs[[cond]] <- cv.nnlvg(X, rhos = logseq(0.01, 0.5, 5), 
                               gammas = logseq(0.01, 0.5, 5))
  nnlvgs[[cond]] <- nnlvg(cov(X), cvnnlvgs[[cond]]$rho, 
                          cvnnlvgs[[cond]]$gamma)
  
  # rcLVGLASSO
  cvrclvgs[[cond]] <- cv.rclvg(X, rhos = logseq(0.01, 0.5, 5), rs = 2:10)
  rclvgs[[cond]] <- rclvg(cov(X), cvrclvgs[[cond]]$rho, 
                          cvrclvgs[[cond]]$r)
  
  # tGLASSO
  tgs[[cond]] <- ebic.tg(X)
  
  # Save
  save(cvslis, slis, 
       cvnnlvgs, nnlvgs,
       cvrclvgs, rclvgs,
       tgs, file = "ests.rda") 
}

### Subject-wise ###
set.seed(123)
slis_all <- vector("list", length(conds))
names(slis_all) <- conds
# Run estimation
for(cond in conds){
  subj_all <- vector("list", 16)
  for(i in 1:16){
    print(paste0("GROUP: ", cond, ", SUBJ: ", i))
    X <- Xs[[cond]][((i-1)*100+1):(i*100),]
    subj_all[[i]] <- slice(cov(X), cvslis[[cond]]$rho,
                           cvslis[[cond]]$r)
  }
  slis_all[[cond]] <- subj_all
}
save(slis_all, file = "subj_ests.rda") # Save