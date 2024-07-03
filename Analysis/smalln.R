library(MASS)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Matrix)
library(scales)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Get file names
plats <- 3:6; ns <- seq(75, 375, 75)
models <- c("SLICE", "nnLVGLASSO", "tGLASSO", "rcLVGLASSO")
allsims <- data.frame()
for(model in models){
  files <- list.files(paste0("../Simulations/", model))
  files <- files[!grepl("n10000", files)]
  for(plat in plats){
    for(n in ns){
      print(paste0("model:", model, ", plat:", plat, ", n:", n))
      currfiles <- files
      currfiles <- currfiles[grepl(paste0("plat", plat), files) & 
                             grepl(paste0("n", n), files)]
      for(j in 1:length(currfiles)){
        load(paste0("../Simulations/", model, "/", currfiles[j]))
        simtype <- strsplit(files[j], "_")[[1]][2]
        
        iters <- which(!sapply(S_stars, is.null))
        F1 <- sin_theta <- frob_norm <- spec_norm <- ARI <- rep(NA, length(iters))
        i <- 1
        for(k in iters){
          # Sparse metrics
          F1[i] <- F1score(S_stars[[k]], S_hats[[k]])
          
          # Low rank metrics
          if(model != "tGLASSO"){
            eigL_hat <- eigen(L_hats[[k]])
            eigL_star <- eigen(L_stars[[k]])
            
            # Sin theta
            sin_theta[i] <- sintheta(eigL_star$vectors[,1], eigL_hat$vectors[,1])
            
            # Frobenius norm
            frob_norm[i] <- norm(L_stars[[k]] - L_hats[[k]])
            
            # The spectral norm
            eigDiff <- eigen(L_stars[[k]] - L_hats[[k]])
            spec_norm[i] <- eigDiff$values[1]
            
            # Do kmeans clustering
            z_hat <- kmeans(eigL_hat$vectors[,1:plat], max(z_stars[[k]]), iter.max = 100, nstart = 1000)$cluster
            ARI[i] <- ari(z_stars[[k]], z_hat)
          }
          i <- i + 1
        }
        # Save results
        allsims <- rbind(allsims, cbind(model, plat, n, iters, F1, sin_theta, frob_norm, spec_norm, ARI)) # Add to df
        save(allsims, file = "lowNdf.rda")
      }
    }
  }
}

# Get summary values
rm(list = setdiff(ls(), "allsims"))
allsims <- allsims %>%
  mutate(plat = as.numeric(plat),
         n = as.numeric(n),
         F1 = as.numeric(F1),
         sin_theta = as.numeric(sin_theta),
         frob_norm = as.numeric(frob_norm),
         spec_norm = as.numeric(spec_norm),
         ARI = as.numeric(ARI))
save(allsims, file = "lowNdf.rda")
load("lowNdf.rda")

# Summarize
df <- allsims %>%
  group_by(model, plat, n) %>%
  summarise(across(c(F1, sin_theta, frob_norm, spec_norm, ARI), mean, na.rm = TRUE))
write.table(df, file = "lowNmean.txt")
sd_df <- allsims %>%
  group_by(model, plat, n) %>%
  summarise(across(c(F1, sin_theta, frob_norm, spec_norm, ARI), sd, na.rm = TRUE))
write.table(sd_df, file = "lowNsd.txt")

# Function to plot and save line graphs
plot_and_save_lines <- function(df, value){
  subdf <- df[,colnames(df) %in% c("model", "plat", "n", value)]
  
  # Change names to shortened version
  subdf$model[subdf$model == "nnLVGLASSO"] = "nnLVG"; subdf$model[subdf$model == "rcLVGLASSO"] = "rcLVG"
  
  # Calculate the mean for each group
  colnames(subdf)[4] <- "value"
  
  # Remove rows for tGLASSO if not calculating sparse
  if(value != "F1"){
    subdf <- subdf[subdf$model != "tGLASSO",]
  }
  
  # Generate ylabs
  if(value == "F1"){
    ylab <- "F1 Score"
    legpos <-  c(0.27, 0.865)
  } else if(value == "sin_theta"){
    ylab <- expression("sin" * Theta(hat(u[1]), u[1]))
    legpos <-  c(0.2, 0.15)
  } else if(value == "spec_norm"){
    ylab <- expression("||" * (hat(L) - L^"*") * "||"[2])
    legpos <-  c(0.2, 0.15)
  } else if(value == "ARI"){
    ylab <- value
    legpos <-  c(0.2, 0.866)
  }
  
  line_types <- c("3" = "solid", "4" = "dashed", "5" = "dotted", "6" = "dotdash")  # Adjust as per your actual plat values
  model_cols <- c("SLICE" = "#0C6291", "nnLVG" = "#A63446", "rcLVG" = "#34A646", "tGLASSO" = "#A634A6")
  p <- ggplot(subdf, aes(x = n, y = value, group = interaction(model, plat), color = model, linetype = as.factor(plat))) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = as.factor(plat)), size = 3.5) +
    labs(x = "n", y = ylab, color = "Model", linetype = "Rank of Latent", shape = "Rank of Latent") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 18, margin = margin(t = 15, r = 0, b = 0, l = 0)), # Adjusted for x-axis title
      axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 5, b = 0, l = 0)), # Adjusted for
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.position = legpos) +
    scale_linetype_manual(values = line_types) +
    scale_color_manual(values = model_cols) +
    scale_x_continuous(breaks = c(75, 150, 225, 300, 375)) +
    guides(
      color = guide_legend(nrow = 1, order = 2),
      linetype = guide_legend(nrow = 1, order = 1),
      shape = guide_legend(nrow = 1, order = 1),
    )
  print(p)
  ggsave(paste0(value, ".png"), plot = p, width = 10, height = 7, dpi = 300)
}

# List of values to loop through
vals <- c("F1", "sin_theta", "spec_norm", "ARI")
for(i in 1:length(vals)){
  plot_and_save_lines(df, vals[i])
}
