library(ggplot2)
library(reshape2)
library(pracma)
load("svmmat.rda")
load("svmpreds.rda")

# Heatmap of confusion matrix
y_pred <- factor(svmpreds$y_pred); y_test <- factor(svmpreds$y_test)
levels(y_pred) <- levels(y_test) <- c("Fam", "Scr", "Unf")

# Generate the confusion matrix
conf_mat <- table(y_test, y_pred)

# Convert the confusion matrix to a data frame
conf_mat_df <- as.data.frame(as.table(conf_mat))
custom_colors <- colorRampPalette(c("white", "#0C6291"))(100)
p <- ggplot(conf_mat_df, aes(x = y_test, y = y_pred, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = custom_colors, limits = c(0, 100)) +
  geom_text(aes(label = Freq), color = "black", size = 20) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 35),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 35),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40), 
    panel.grid = element_blank(),
    legend.position = "none",
    aspect.ratio = 1,
  ) +
  labs(title = "", x = "Predicted", y = "True")
p
ggsave(paste0("confmat.png"), plot = p, width = 8, height = 8, dpi = 300)