library(MASS)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Generate a covariance matrix
set.seed(123)
invSigma = matrix(rnorm(200), 10, 20)
invSigma = invSigma %*% t(invSigma)
Sigma = solve(invSigma)

# Check likelihood
mle = logL(Sigma, solve(Sigma))

# Look at different samples
ns = c(100, 500, 1000, 5000, 10000, 50000, 100000, 500000)

# Loop through different sample sizes
delta_logL = matrix(NA, 100, length(ns))
delta_inv = matrix(NA, 100, length(ns))
for(i in 1:length(ns)){
  for(j in 1:100){
    covmat = cov(mvrnorm(ns[i], rep(0, 10), Sigma))
    delta_logL[j,i] = mle - logL(Sigma, covmat)
    delta_inv[j,i] = norm(solve(covmat) - invSigma, "F")
  }
}

# Plot the means
library(ggplot2)
library(dplyr)
library(tidyr)

# Create a data frame from your matrices
data <- data.frame(ns = rep(ns, each = 100), 
                   delta_logL = c(delta_logL), 
                   delta_inv = c(delta_inv))

# Calculate means and standard deviations
data_summary <- data %>%
  group_by(ns) %>%
  summarize(mean_logL = mean(delta_logL), sd_logL = sd(delta_logL),
            mean_inv = mean(delta_inv), sd_inv = sd(delta_inv)) %>%
  ungroup()

# Pivot longer for ggplot
data_long <- pivot_longer(data_summary, 
                          c(mean_logL, sd_logL, mean_inv, sd_inv),
                          names_to = "metric", values_to = "value")

# Separate metric and statistic into separate columns
data_long <- data_long %>%
  separate(metric, into = c("statistic", "metric"), sep = "_")

# Filter data for means and standard deviations separately
mean_data <- filter(data_long, statistic == "mean")
sd_data <- filter(data_long, statistic == "sd")

# Merge mean and sd data to include sd as a column in mean_data
plot_data <- merge(mean_data, sd_data, by = c("ns", "metric"))

# Rename columns for clarity
colnames(plot_data) <- c("ns", "metric", "statistic_mean", "mean_value", "statistic_sd", "sd_value")

# Plot
ggplot(plot_data, aes(x = ns, y = mean_value, color = metric)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, fill = metric), alpha = 0.2) +
  scale_x_continuous(name = "Sample Size") +
  scale_y_continuous(name = "Value") +
  labs(color = "Metric", fill = "Metric") +
  theme_minimal()
