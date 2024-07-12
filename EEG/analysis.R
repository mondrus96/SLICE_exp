sapply((paste0("../Core/", list.files("../Core/"))), source)
load("ests.rda")

# Look through slis
heatmap(slis[[1]]$L, Rowv = NA, Colv = NA)
heatmap(slis[[2]]$L, Rowv = NA, Colv = NA)
heatmap(slis[[3]]$L, Rowv = NA, Colv = NA)

# Look at 
sapply(cvslis, "[[", "r")
