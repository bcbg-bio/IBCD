library(caret)
library(readr)
library(pcalg)
library(parallel)

# =========================================================
# helper: fit GES on train split
# =========================================================
fit_ges_on_train <- function(Y_train, y_train, lam_base) {
  n_train <- nrow(Y_train)
  lam <- lam_base * log(n_train)
  
  score_obj_train <- new(
    "GaussL0penObsScore",
    data = as.matrix(Y_train),
    lambda = lam
  )
  
  t0 <- Sys.time()
  ges_res <- ges(score_obj_train)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  list(
    ges_res = ges_res,
    lambda = lam,
    elapsed = elapsed
  )
}

# =========================================================
# helper: held-out score on test split
# IMPORTANT:
# fit graph on train, then evaluate that graph on held-out test
# using a new observational score object built from test data
# =========================================================
heldout_global_score_ges <- function(Y_test, y_test, dag_repr, lam_base) {
  n_test <- nrow(Y_test)
  lam_test <- lam_base * log(n_test)
  
  score_obj_test <- new(
    "GaussL0penObsScore",
    data = as.matrix(Y_test),
    lambda = lam_test
  )
  
  score_obj_test$global.score(dag_repr)
}

# =========================================================
# load observational dataset
# =========================================================
base_dir <- "/path/to/data"

out_dir <- file.path(base_dir, "500d", "sf", "42")

X_fn <- file.path(out_dir, "Y_matrix.csv")
y_fn <- file.path(out_dir, "targets.txt")

Y <- read.csv(X_fn, check.names = FALSE)
target_strings <- readLines(y_fn)
stopifnot(nrow(Y) == length(target_strings))

# =========================================================
# make tuning folds from all rows
# seed 456 = tuning split
# =========================================================
set.seed(456)

folds_all <- caret::createFolds(
  seq_len(nrow(Y)),
  k = 5,
  returnTrain = FALSE
)

# =========================================================
# lambda grid: 10 log-spaced values from 0.5 to 64
# =========================================================
lam_grid <- exp(seq(log(0.5), log(64), length.out = 10))

all_results <- list()

for (g in seq_along(lam_grid)) {
  lam_base <- lam_grid[g]
  
  cat(sprintf(
    "\n================ LAMBDA %d/%d | lam_base=%.6f ================\n",
    g, length(lam_grid), lam_base
  ))
  flush.console()
  
  fold_scores <- numeric(1)
  fold_times <- numeric(1)
  
  # quick version: only one held-out fold
  for (fold in 1:5) {
    te_idx <- folds_all[[fold]]
    tr_idx <- setdiff(seq_len(nrow(Y)), te_idx)
    
    Y_tr <- Y[tr_idx, , drop = FALSE]
    y_tr <- target_strings[tr_idx]
    
    Y_te <- Y[te_idx, , drop = FALSE]
    y_te <- target_strings[te_idx]
    
    cat(sprintf(
      "[fold %d/5] train=%d test=%d | lam_base=%.6f\n",
      fold, nrow(Y_tr), nrow(Y_te), lam_base
    ))
    flush.console()
    
    fit <- fit_ges_on_train(Y_tr, y_tr, lam_base)
    
    score_val <- heldout_global_score_ges(
      Y_test = Y_te,
      y_test = y_te,
      dag_repr = fit$ges_res$repr,
      lam_base = lam_base
    )
    
    fold_scores[fold] <- score_val
    fold_times[fold] <- fit$elapsed
    
    cat(sprintf(
      "   done | heldout_global_score_ges=%.6f | ges_time=%.1fs\n",
      score_val, fit$elapsed
    ))
    flush.console()
  }
  
  mean_score <- mean(fold_scores)
  sd_score <- NA
  
  cat(sprintf(
    ">>> lam_base=%.6f | mean heldout score=%.6f\n",
    lam_base, mean_score
  ))
  flush.console()
  
  all_results[[g]] <- data.frame(
    lam_base = lam_base,
    fold1_score = fold_scores[1],
    mean_heldout_score = mean_score,
    sd_heldout_score = sd_score,
    mean_ges_time_sec = mean(fold_times)
  )
}

results_df <- do.call(rbind, all_results)
results_df <- results_df[order(-results_df$mean_heldout_score), ]

cat("\n================ FINAL TUNING RESULTS ================\n")
print(results_df)

best_lam_base <- results_df$lam_base[1]
cat(sprintf("\nBest global lam_base = %.6f\n", best_lam_base))