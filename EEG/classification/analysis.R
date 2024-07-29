library(ggplot2)
library(reshape2)
library(pracma)
load("svmmats.rda")
load("svmouts.rda")

# Hypothesis testing
prop.test(alternative = "greater")

# Heatmap of confusion matrix
for(model in names(svmouts)){
  y_pred <- factor(svmouts[[model]]$y_pred)
  y_test <- factor(svmouts[[model]]$y_test)
  
  # Generate the confusion matrix
  conf_mat <- table(y_test, y_pred)
  
  # Plot the heatmap with the specified customizations
  if(model == "slice"){
    model_nme <- "SLICE"
  } else if(model == "rclvg"){
    model_nme <- "rcLVG"
  } else if(model == "nnlvg"){
    model_nme <- "nnLVG"
  }
  # Convert the confusion matrix to a data frame
  conf_mat_df <- as.data.frame(as.table(conf_mat))
  custom_colors <- colorRampPalette(c("white", "#0C6291"))(100)
  p <- ggplot(conf_mat_df, aes(x = y_test, y = y_pred, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = custom_colors, limits = c(0, 100)) +
    geom_text(aes(label = Freq), color = "black", size = 15) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
      axis.text.y = element_text(angle = 90, hjust = 0.5, size = 20),
      panel.grid = element_blank(),
      legend.position = "none",
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, size = 30)
    ) +
    labs(title = model_nme, x = "", y = "")
  ggsave(paste0(model, "_confmat.png"), plot = p, width = 8, height = 8, dpi = 300)
}