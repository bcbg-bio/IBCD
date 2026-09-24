library(pcalg)
library(readr)
library(SID)
library(dplyr)   
library(purrr)   

calc_metrics <- function(X, X_true, eps = 1e-8) {
  D <- ncol(X)
  X <- X[!diag(D)]
  X_true <- X_true[!diag(D)]
  
  X[abs(X) < eps] <- 0
  X_true[abs(X_true) < eps] <- 0
  
  rmse <- sqrt(mean((X - X_true)^2))
  mae <- mean(abs(X - X_true))
  
  sign_X <- sign(X)
  sign_Xt <- sign(X_true)
  TN <- sum(!(abs(sign_X) | abs(sign_Xt)))
  TS <- sum(sign_X == sign_Xt) - TN
  FN <- sum((1 - abs(sign_X)) & abs(sign_Xt))
  FS <- sum(sign_X != sign_Xt) - FN
  
  shd <- sum(sign_X != sign_Xt)
  N_pos <- sum(X_true > 0)
  N_neg <- sum(X_true < 0)
  N_zero <- sum(X_true == 0)
  
  weights <- sign_Xt
  weights[sign_Xt > 0] <- 1 / N_pos
  weights[sign_Xt < 0] <- 1 / N_neg
  weights[sign_Xt == 0] <- 1 / N_zero
  weight_acc <- sum((sign_X == sign_Xt) * weights) / sum(weights)
  
  precision <- TS / (TS + FS)
  recall <- TS / (TS + FN)
  f1 <- 2 * precision * recall / (precision + recall)
  
  return(list(precision = precision, recall = recall, F1 = f1,
              rmse = rmse, mae = mae, shd = shd, weight_acc = weight_acc,
              TP = TS, FP = FS, TN = TN, FN = FN))
}


# ---- settings ----
dir_path <- "/data/path"   
#files <- list.files(dir_path, pattern = "G_proper_*\\.csv", full.names = TRUE) #posterior_mean_sf_essential_fold1.csv
files <- list.files(dir_path, pattern = "G_.*\\.csv", full.names = TRUE)
eps_val <- 0.075

# ---- load all matrices ----
G_list <- lapply(files, function(f) {
  as.matrix(read.csv(f, row.names = 1))
})

names(G_list) <- basename(files)

# ---- metric function (your existing one assumed) ----
# calc_metrics(est, true, eps)

# ---- pairwise comparisons ----
results <- list()
k <- 1

for (i in 1:(length(G_list)-1)) {
  for (j in (i+1):length(G_list)) {
    
    G1 <- G_list[[i]]
    G2 <- G_list[[j]]
    
    m <- calc_metrics(abs(G1), abs(G2), eps = eps_val)
    
    results[[k]] <- data.frame(
      pair = paste(names(G_list)[i], names(G_list)[j], sep = " vs "),
      SHD = m$shd,
      F1  = m$F1
    )
    
    k <- k + 1
  }
}

results_df <- do.call(rbind, results)

print(results_df)

# ---- summary ----
mean_shd <- mean(results_df$SHD)
sd_shd   <- sd(results_df$SHD)

mean_f1 <- mean(results_df$F1)
sd_f1   <- sd(results_df$F1)

cat(sprintf("\nSHD: %.3f ± %.3f\n", mean_shd, sd_shd))
cat(sprintf("F1 : %.3f ± %.3f\n", mean_f1, sd_f1))



