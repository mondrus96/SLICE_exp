library(ggplot2)
library(dplyr)
library(Matrix)
library(gridExtra)

# For plotting example Ls recovered by different methods
sapply((paste0("../Core/", list.files("../Core/"))), source)

### Plotting Examples ###
set.seed(123)
pobs <- 150 # Number of observed variables for S
n <- 10000 # Number of observations

Lnames <- c("cres", "spir")
models <- c("True", "rcLVG", "nnLVG", "SLICE")
Lout <- list(Lcres(pobs, 0.1), Lspir(pobs, 0.02))
S_star <- Smat(pobs, 2, 1.5); S_star[S_star < 0.01] <- 0
for(i in 1:length(Lout)){
  set.seed(123)
  L_star <- Lout[[i]]$L
  Sigma_star <- solve(S_star + L_star)
  X <- mvrnorm(n, rep(0, pobs), Sigma_star)
  Sigma <- cov(X)
  
  cvsli <- cv.slice(X); cvnnlvg <- cv.nnlvg(X); cvrclvg <- cv.rclvg(X)
  
  sliout <- slice(Sigma, cvsli$lambda, cvsli$r)
  nnlvgout <- nnlvg(Sigma, cvnnlvg$lambda, cvnnlvg$gamma)
  rclvgout <- rclvg(Sigma, cvrclvg$lambda, cvrclvg$r)
  
  Lall <- list(eigen(L_star)$vectors[,1:2], 
               eigen(rclvgout$L)$vectors[,1:2],
               eigen(nnlvgout$L)$vectors[,1:2],
               eigen(sliout$L)$vectors[,1:2])
  plotlist <- vector("list", 4)
  for(j in 1:length(models)){
    currplot <- ggplot(as.data.frame(Lall[[j]]), aes(x = V1, y = V2)) +
      geom_point() + theme_minimal() + geom_point(colour = "black", alpha = 0.5) +
      labs(title = models[j], x = NULL, y = NULL) +  # Set x and y labels to NULL
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_blank(),  # Remove x-axis tick labels
        axis.text.y = element_blank(),  # Remove y-axis tick labels
        axis.title.x = element_blank(),  # Optionally, remove x-axis title
        axis.title.y = element_blank(),  # Optionally, remove y-axis title
        axis.ticks = element_blank()     # Optionally, remove axis ticks
      ) +
      coord_fixed(ratio = 1) +
      xlim(c(min(Lall[[1]][,1]), max(Lall[[1]][,1]))) + 
      ylim(c(min(Lall[[1]][,2]), max(Lall[[1]][,2])))
    plotlist[[j]] <- currplot
  }
  # Arrange the plots in a grid
  gridplot <- grid.arrange(grobs = plotlist, ncol = 4)
  ggsave(paste0(Lnames[i], ".png"), plot = gridplot, width = 8, height = 3, dpi = 600)
}