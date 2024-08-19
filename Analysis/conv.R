library(MASS)
library(ggplot2)
sapply((paste0("../Core/", list.files("../Core/"))), source)

set.seed(123)

# Simulate data at different values of rank
pobs <- 200 # Number of observed variables for S
rs <- 2:6

# Loop through different rs
deltaL <- deltaS <- deltalogL <- vector("list", length(rs))
for(i in 1:length(rs)){
  print(rs[i])
  
  S_star <- Smat(pobs, 2, 1.5)
  S_star[S_star < 0.01] <- 0
  L_star <- Lrand(pobs, rs[i], 1.5)$L
  Sigma <- solve(S_star + L_star)
  
  sliout <- slice(Sigma, 0.01, rs[i])
  deltaS[[i]] <- sliout$misc$deltaS; deltaL[[i]] <- sliout$misc$deltaL; deltalogL[[i]] <- sliout$misc$deltalogL
}

plot_and_save <- function(data_list, file_name, xlab, ylab) {
  # Transform the list into a long-format data frame
  data_long <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    data.frame(Iteration = seq_along(data_list[[i]]),
               Value = data_list[[i]],
               Line = factor(i))
  }))
  
  # Create the plot
  p <- ggplot(data_long, aes(x = Iteration, y = Value, group = Line, color = Line)) +
    geom_line(linewidth = 1.5) +
    labs(x = xlab, y = ylab, color = "Rank of Latent") +
    theme_bw() +
    theme(
      legend.position = c(0.85, 0.82),
      legend.title = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      legend.text = element_text(size = 25))
  print(p)
  
  # Save the plot as a PNG file
  png(filename = file_name, width = 2000, height = 1400, res = 300)
  print(p)
  dev.off()
}

# Now, call the function for each of your lists
plot_and_save(deltaS, "deltaS_plot.png", "Iterations", expression(Delta * S))
plot_and_save(deltaL, "deltaL_plot.png", "Iterations", expression(Delta * L))
plot_and_save(deltalogL, "deltalogL_plot.png", "Iterations", expression(Delta * log-Liklihood))
