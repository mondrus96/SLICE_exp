library(ggplot2)
library(dplyr)
library(Matrix)

# For plotting example Ls recovered by different methods
sapply((paste0("../Core/", list.files("../Core/"))), source)

### Plotting Examples ###
set.seed(123)
pobs <- 150 # Number of observed variables for S
n <- 10000 # Number of observations

Lnames <- c("cres", "spir")
models <- c("True", "SLICE", "nnLVGLASSO", "rcLVGLASSO")
Lout <- list(Lcres(pobs, 0.1), Lspir(pobs, 0.02))
S_star <- Smat(pobs, 2, 1.5); S_star[S_star < 0.01] <- 0
for(i in 1:length(Lout)){
  L_star <- Lout[[i]]$L
  Sigma_star <- solve(S_star + L_star)
  X <- mvrnorm(n, rep(0, pobs), Sigma_star)
  Sigma <- cov(X)
  
  cvsli <- cv.slice(X); cvnnlvg <- cv.nnlvg(X); cvrclvg <- cv.rclvg(X)
  
  sli <- slice(Sigma, cvsli$lambda, cvsli$r)
  nnlvg <- nnlvg(Sigma, cvnnlvg$lambda, cvnnlvg$gamma)
  rclvg <- rclvg(Sigma, cvrclvg$lambda, cvrclvg$r)
  
  Lall <- list(eigen(L_star)$vectors[,1:2], 
               eigen(sli$L)$vectors[,1:2],
               eigen(nnlvg$L)$vectors[,1:2],
               eigen(rclvg$L)$vectors[,1:2])
  for(j in 1:length(models)){
    currplot <- ggplot(as.data.frame(Lall[[j]]), aes(x = V1, y = V2)) +
      geom_point() + theme_minimal() + geom_point(colour = "blue", alpha = 0.5) +
      labs(title = models[j],
           x = "",
           y = "") +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0(Lnames[i], "_", models[j], ".png"), plot = currplot, width = 8, height = 6, dpi = 600)  
  }
}
