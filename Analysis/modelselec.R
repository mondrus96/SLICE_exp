library(ggplot2)
library(tidyr)
library(dplyr)
sapply((paste0("../Core/", list.files("../Core/"))), source)

### Model Selection ###
# Simulation parameters
pobs <- 150 # Number of observed variables for S
plat <- 7 # Number of latent variables for L
n <- 10000 # Number of observations
simtype <- "rand"
iters <- 100
rhoseq <- logseq(5e-3, 0.05, 10)
rseq <- 2:11

rs <- rhos <- c()
likls <- matrix(0, length(rseq), length(rhoseq))
for(i in 1:iters){
  set.seed(123*i)
  
  # Generate simulation data
  print(paste0("SIM ITER ", i))
  S_star <- Smat(pobs, 2, 1.5)
  S_star[S_star < 0.01] <- 0 # True sparse component
  Lout <- Lrand(pobs, plat, 1.5)
  L_star <- Lout$L # True latent component
  
  Sigma_star <- solve(S_star + L_star) # True Sigma
  X <- mvrnorm(n, rep(0, pobs), Sigma = Sigma_star) # Finite sample data
  Sigma <- cov(X) # Sample Sigma 
  
  cvsli <- cv.slice(X, 3, rhoseq, rseq)
  rs <- c(rs, cvsli$r); rhos <- c(rhos, cvsli$rho)
  likls <- likls + exp(cvsli$cvmat)
}
likls <- likls/iters
# Save
save(rs, rhos, likls, file = "modelselec.rda")
load("modelselec.rda")

# Convert it to long format
likl_long <- as.data.frame(likls) %>%
  tibble::rownames_to_column("r") %>%
  gather(key = "rho", value = "likl", -r) %>%
  mutate(rho = factor(rho, levels = sort(unique(as.numeric(rho)))),
         r = factor(r, levels = sort(as.numeric(unique(r)))))
selec_df <- data.frame(r = factor(rs, levels = sort(unique(as.numeric(rs)))), 
                         rho = factor(rhos, levels = sort(unique(as.numeric(rhos)))))
total_count <- nrow(selec_df)
selec_count <- selec_df %>%
  group_by(r, rho) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / total_count)
likl_long <- merge(likl_long, selec_count, by = c("r", "rho"), all.x = TRUE)

# For formatting
format_rho <- function(x) {
  sprintf("%.4f", as.numeric(x))
}

# Plot
p <- ggplot() +
  geom_tile(data = likl_long, aes(x = rho, y = r, fill = likl), na.rm = TRUE) +
  geom_text(data = subset(likl_long, !is.na(proportion)), 
            aes(x = rho, y = r, label = sprintf("%.2f", proportion)), 
            color = "blue", size = 8, na.rm = TRUE) +
  scale_x_discrete(labels = format_rho) +
  scale_fill_gradient(low = "#636363", high = "white") +  # Set color gradient from black to white
  labs(x = "rho", y = "Rank") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position = "none"  # Remove the legend
  )
print(p)

# Display the plot
ggsave("modelselec.png", p, width = 12, height = 7, units = "in")
