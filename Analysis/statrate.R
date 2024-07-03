library(MASS)
library(ggplot2)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Define simulation parameters
set.seed(321)
p <- 100
L_star <- Lrand(p, 2, 0.4)$L
S_star <- Smat(p, 1, 0.8)
sort(unique(as.vector(S_star)), TRUE)
S_star[S_star < 0.2] <- 0 # True sparse component
s <- S_star
diag(s) <- 0
s <- sum(s != 0) # Number of non-zero edges
print(s)
Sigma_star <- solve(S_star + L_star)

# Loop through different values of n
ns <- seq(200, 1000, 5)
L_theo_r <- S_theo_r <- Ldist <- Sdist <- c()
lambda <- 0.2
for(i in seq_along(ns)){
  print(i)
  n <- ns[i] # Current sample size
  
  L_theo_r <- c(L_theo_r, sqrt(log(p)/n)) # Theoretical rate of L
  S_theo_r <- c(S_theo_r, sqrt(s*log(p)/n)) # Theoretical rate of S
  
  X <- mvrnorm(n, rep(0, p), Sigma_star) # Finite sample simulation
  Sigma <- cov(X)
  
  sliout <- slice(Sigma, lambda*sqrt(log(p)/n), 2, tol = 1e-3, maxiter = 1000) # Estimates
  
  Ldist <- c(Ldist, norm(sliout$L - L_star, type = "I")) # Actual distance in L
  S_off <- sliout$S - S_star
  diag(S_off) <- 0
  Sdist <- c(Sdist, norm(S_off, type = "F")) # Actual distance in S
}
L_obs_r <- Ldist/mean(Ldist/L_theo_r)
S_obs_r <- Sdist/mean(Sdist/S_theo_r)

# Assuming Lbound and adjusted_observed are vectors of the same length
data <- data.frame(
  n = ns,
  Rate = c(L_theo_r, L_obs_r, S_theo_r, S_obs_r),
  Type = rep(c("L Theoretical", "L Estimated", "S Theoretical", "S Estimated"), each = length(L_theo_r))
)

# Plotting only S
data_s <- data[grep("^S", data$Type),]
png("Srate.png", width = 8, height = 5, units = "in", res = 600, pointsize = 10)
ggplot(data_s, aes(x = n, y = Rate, color = Type, group = Type, shape = Type)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c("S Theoretical" = "black", "S Estimated" = "#0C6291")) +
  labs(title = expression("||" * hat(S)[off] - S[off]^"*" * "||"[F]), x = "n", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")
dev.off()

# Plotting only L
data_l <- data[grep("^L", data$Type),]
png("Lrate.png", width = 8, height = 5, units = "in", res = 600, pointsize = 10)
ggplot(data_l, aes(x = n, y = Rate, color = Type, group = Type, shape = Type)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c("L Theoretical" = "black", "L Estimated" = "#0C6291")) +
  labs(title = expression("||" * hat(L) - L^"*" * "||"[infinity]), x = "n", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none")
dev.off()
