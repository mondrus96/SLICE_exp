library(MASS)
library(ggplot2)
library(dplyr)
sapply((paste0("../Core/", list.files("../Core/"))), source)

# Get file names
files <- list.files("../Simulations", "n10000")

# Loop through files
allsims <- c()
for(i in 1:length(files)){
  simtype <- strsplit(strsplit(files[i], "_")[[1]], "sim")[[1]][2]
  df <- cbind(simtype, read.table(paste0("../Simulations/", files[i]), header = TRUE))
  df <- df[!colnames(df) %in% c("lvg_nmi", "sli_nmi", "lvg_ari", "sli_ari")]
  allsims <- rbind(allsims, df)
}

# Get summary values
dfsummary <- allsims %>%
  group_by(simtype) %>%
  summarize(across(c(lvg_sin, sli_sin, lvg_fnorm, sli_fnorm), 
                   list(mean = ~ mean(., na.rm = TRUE), 
                        sd = ~ sd(., na.rm = TRUE)), 
                   .names = "{.col}_{.fn}"))

# Example plots
set.seed(123)
pobs <- 150 # Number of observed variables for S
n <- 10000 # Number of observations

S_star <- Smat(pobs, 2, 1.5)
Lout <- list(Lspir(pobs, 0.05), Lcirc(pobs, 0.05), Lcres(pobs, 0.1))
Lnames <- c("spir", "circ", "cres")
models = c("True", "slice", "lvglasso")
for(i in 1:length(Lout)){
  L_star <- Lout[[i]]$L
  Sigma_star <- solve(S_star + L_star)
  X <- mvrnorm(n, rep(0, pobs), Sigma_star)
  Sigma <- cov(X)
  
  cvsli <- cv.slice(X)
  cvlvg <- cv.lvg(X)
  
  sliout <- slice(Sigma, cvsli$lambda, cvsli$r)
  lvgout <- lvglasso(Sigma, cvlvg$lambda, cvlvg$gamma) 
  
  Lall <- list(eigen(L_star)$vectors[,1:2], 
               eigen(sliout$L)$vectors[,1:2],
               eigen(lvgout$L)$vectors[,1:2])
  for(j in 1:length(models)){
    currplot <- ggplot(as.data.frame(Lall[[j]]), aes(x = V1, y = V2)) +
      geom_point() + theme_bw() + geom_point(colour = "blue", alpha = 0.5) +
      labs(title = models[j],
           x = "",
           y = "") +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0(Lnames[i], "_", models[j], ".png"), plot = currplot, width = 8, height = 6, dpi = 600)  
  }
}
