library(MASS)
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Get file names
files <- list.files("../Simulations", "_n")
files <- files[!grepl("n10000", files)]

# Check which files are missing
rank <- c(3, 4, 5, 6)
n <- c(75, 150, 225, 300, 375)
allsims <- c()
for(plat in rank){
  for(num in n){
    for(i in 1:4){
      start <- (1 + 25 * (i - 1)); end <- 25 * i
      allsims <- c(allsims, paste0("simrand_plat", plat, "_n", 
                                   num, "_iters", start, "to", end, ".txt"))
    }
  }
}

# Didn't run
print(allsims[!(allsims %in% files)])

# Loop through files
df <- c()
incomp <- c()
for(i in 1:length(files)){
  sim <- read.table(paste0("../Simulations/", files[i]), header = TRUE)
  if(nrow(sim) < 25){
    incomp <- c(incomp, files[i])
  }
  df <- rbind(df, sim)
}
# Started, but didn't complete
print(incomp)

# Summarize
df <- df %>%
  group_by(plat, n) %>%
  summarise(across(c(lvg_nmi, lvg_ari, lvg_sin, lvg_fnorm, 
                     sli_nmi, sli_ari, sli_sin, sli_fnorm), mean, na.rm = TRUE))

# Function to plot and save heatmaps
plot_and_save_heatmap <- function(data, value, filename_prefix, scale_limits){
  data <- data[,colnames(data) %in% c("plat", "n", value)]
  colnames(data)[3] <- "value"
  p <- ggplot(data, aes(plat, n, fill = value)) +
    geom_tile() + scale_fill_viridis(limits = scale_limits, alpha = 0.8) + 
    geom_text(aes(label = round(value, 2)), vjust = 0.5, hjust = 0.5) +
    labs(x = "plat", y = "n", title = paste0(value)) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),                       
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_blank())
  
  ggsave(paste0(value, ".png"), plot = p, width = 8, height = 6, dpi = 300)
}

# Function to plot and save line grahps
plot_and_save_lines <- function(df, value){
  subdf <- df %>%
    select(plat, n, contains(value))
  
  # Reshape data to long format
  longsubdf <- subdf %>%
    pivot_longer(cols = starts_with(paste0("lvg_", value)) | starts_with(paste0("sli_", value)), 
                 names_to = "method", values_to = value) %>%
    mutate(method = gsub(paste0("_", value), "", method))
  
  # Calculate the mean for each group
  colnames(longsubdf)[4] <- "value"
  meansubdf <- longsubdf %>%
    group_by(plat, n, method) %>%
    summarise(mean_val = mean(value, na.rm = TRUE))
  
  # Create the plot
  p <- ggplot(meansubdf, aes(x = n, y = mean_val, group = interaction(method, plat), color = as.factor(plat), linetype = method)) +
    geom_line() +
    geom_point(aes(shape = method)) +
    scale_shape_manual(values = c("lvg" = 15, "sli" = 17)) +  # 15: filled square, 17: filled triangle
    scale_linetype_manual(values = c("lvg" = "solid", "sli" = "dashed")) +
    labs(x = "n", y = vals[i], color = "Plat", linetype = "Method", shape = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(value, ".png"), plot = p, width = 8, height = 6, dpi = 300)
}

# List of values to loop through
values <- c("nmi", "ari", "sin", "fnorm")
limits_list <- list(c(0, 1), c(0, 1), c(0, 1), c(0, 200))

# Loop through each method and each value
for(method in c("lvg", "sli")) {
  for(i in seq_along(values)) {
    value <- paste(method, values[i], sep = "_")
    plot_and_save_heatmap(df, value, method, limits_list[[i]])
  }
}

vals <- c("nmi", "ari", "sin", "fnorm")
for(i in 1:length(vals)){
  plot_and_save_lines(df, vals[i])
}