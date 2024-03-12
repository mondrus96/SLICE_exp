library(MASS)
library(ggplot2)
library(reshape2)
sapply((paste0("../Core/", list.files("../Core/"))), source)

### Plot and save TRUE
pobs <- 100 # Number of observed variables for S
plat <- 4 # Number of latent variables for L
n <- 500 # Number of observations
simtype <- "rand"
iters <- 1:100

S_star <- Smat(pobs, 2, 1.5)
S_star[S_star < 0.01] <- 0 # True sparse component
S_star <- 1*(S_star != 0)

Lout <- Lrand(pobs, plat, 1.5)
L_star <- Lout$L; z_star <- Lout$z # True latent component; True cluster labels

# Plot the heatmap in black and white
png("S_star.png", width = 600, height = 600)
heatmap(x = S_star,
        Rowv = NA,  # Disable row clustering
        Colv = NA,  # Disable column clustering
        symm = TRUE, # Indicates the matrix is symmetric
        scale = "none", # Don't scale the data
        labRow = "", # Remove x-axis label
        labCol = "", # Remove y-axis label
        col = colorRampPalette(c("white", "black"))(256)) # Black & White color scheme
dev.off()
png("L_star.png", width = 600, height = 600)
heatmap(x = L_star,
        Rowv = NA,  # Disable row clustering
        Colv = NA,  # Disable column clustering
        symm = TRUE, # Indicates the matrix is symmetric
        scale = "none", # Don't scale the data
        labRow = "", # Remove x-axis label
        labCol = "", # Remove y-axis label
        col = colorRampPalette(c("white", "black"))(256)) # Black & White color scheme
dev.off()

### Plot and save estimated
# Load data
load("../Simulations/diffSest.rda")
mainlist$glasso[[1]]$S
for(model in names(mainlist)){
  # Calculate average S and L
  avgS <- avgL <- matrix(0, 100, 100)
  for(i in 1:length(mainlist[[model]])){
    S <- mainlist[[model]][[i]]$S
    avgS <- avgS + 1*(S != 0)
    avgL <- avgL + mainlist[[model]][[i]]$L
  }
  avgS <- avgS/100; avgL <- avgL/100
  
  # Plot the heatmap in black and white
  png(paste0("S", model, ".png"), width = 600, height = 600)
  heatmap(x = avgS,
          Rowv = NA,  # Disable row clustering
          Colv = NA,  # Disable column clustering
          symm = TRUE, # Indicates the matrix is symmetric
          scale = "none", # Don't scale the data
          labRow = "", # Remove x-axis label
          labCol = "", # Remove y-axis label
          col = colorRampPalette(c("white", "black"))(256)) # Black & White color scheme
  dev.off()
  png(paste0("L", model, ".png"), width = 600, height = 600)
  heatmap(x = avgL,
          Rowv = NA,  # Disable row clustering
          Colv = NA,  # Disable column clustering
          symm = TRUE, # Indicates the matrix is symmetric
          scale = "none", # Don't scale the data
          labRow = "", # Remove x-axis label
          labCol = "", # Remove y-axis label
          col = colorRampPalette(c("white", "black"))(256)) # Black & White color scheme
  dev.off()
}
