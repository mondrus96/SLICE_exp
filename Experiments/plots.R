library(Rtsne)
library(fabisearch)
library(ggplot2)

# Load data
Xall <- read.table("visit2.txt")
AAL <- read.table("AAL.txt")

# Get first subjets data
X <- Xall[1:190,]

# Get estimate
cvout <- cv.slice(X, lambdas = seq(0.1, 0.5, 0.1)) # SLICE
out <- slice(cov(X), cvout$lambda, cvout$r)

# For Plotting
t_sne1 <- Rtsne(data, perplexity = 30, theta = 0.5, dims = 2)
