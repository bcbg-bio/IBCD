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


parse_seeds <- function(s) {
  s <- trimws(s)
  if (grepl(":", s)) {
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    return(seq(as.integer(parts[1]), as.integer(parts[2])))
  } else {
    return(as.integer(strsplit(s, ",", fixed = TRUE)[[1]]))
  }
}

get_arg <- function(args, key, default = NULL) {
  idx <- match(key, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

args <- commandArgs(trailingOnly = TRUE)

g_true <- get_arg(args, "--g_true")
g_est  <- get_arg(args, "--g_est")
method_name  <- get_arg(args, "--method", "unknown")
out_csv       <- get_arg(args, "--out_csv", "summary.csv")
seeds_str     <- get_arg(args, "--seeds", "42:51")
eps_val       <- as.numeric(get_arg(args, "--eps", "0.075"))

if (is.null(g_true) || is.null(g_est)) {
  stop("Need --g_true and --g_est (use {seed} placeholder).")
}

seeds <- parse_seeds(seeds_str)
rows <- vector("list", length(seeds))

for (i in seq_along(seeds)) {
  seed <- seeds[i]
  
  true_path <- gsub("\\{seed\\}", seed, g_true)
  est_path  <- gsub("\\{seed\\}", seed, g_est)
  
  G_true_mat <- as.matrix(read.csv(true_path))
  G_est_mat  <- as.matrix(read.csv(est_path))
  
  m <- calc_metrics(abs(G_est_mat), abs(G_true_mat), eps = eps_val)
  
  A_true <- 1 * (abs(G_true_mat) > eps_val)
  A_est  <- 1 * (abs(G_est_mat)  > eps_val)
  sid_val <- structIntervDist(A_true, A_est)$sid
  
  rows[[i]] <- data.frame(
    seed = seed,
    precision = m$precision,
    recall = m$recall,
    F1 = m$F1,
    rmse = m$rmse,
    mae = m$mae,
    shd = m$shd,
    weight_acc = m$weight_acc,
    SID = sid_val
  )
}

df <- bind_rows(rows)

summary_df <- df %>%
  summarise(
    across(
      c(precision, recall, F1, rmse, mae, shd, weight_acc, SID),
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    )
  )

summary_df <- summary_df %>%
  mutate(method = method_name) %>%
  select(method, everything())


if (file.exists(out_csv)) {
  write.table(
    summary_df,
    out_csv,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
} else {
  write.csv(summary_df, out_csv, row.names = FALSE)
}

cat("Updated summary:", out_csv, "\n")

