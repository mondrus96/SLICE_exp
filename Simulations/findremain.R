# All combos
models <- c("SLICE", "nnLVGLASSO", "tGLASSO", "rcLVGLASSO")
rs <- c(3, 4, 5, 6)
ns <- c(75, 150, 225, 300, 375)
batch <- 1:4

# Loop through all folders and read each sim
df <- c()
nostart <- c()
nofinish <- c()
for(model in models){
  files <- list.files(model)
  for(r in rs){
    for(n in ns){
      for(iter in batch){
        # Check if exists
        start <- (1 + (25*(iter - 1))); end <- (25*iter)
        filename <- paste0(model, "_simrand_plat", r, "_n", n, "_iters", start, "to", end, ".rda")
        if(!filename %in% files){
          nostart <- c(nostart, filename)
          df <- rbind(df, c(model, r, n, iter))
        } else{
          load(paste0(model, "/", filename))
          if(sum(!sapply(S_hats, is.null)) < 25){
            nofinish <- c(nofinish, filename)
            df <- rbind(df, c(model, r, n, iter))
          } 
        }
      }
    }
  }
}
# Save as table
write.table(df, file = "remain.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
