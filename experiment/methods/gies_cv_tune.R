library(caret)
library(readr)
library(pcalg)
library(parallel)


#X_fn <- file.path(out_dir, "Y_matrix.csv")
#y_fn <- file.path(out_dir, "targets.txt")

#Y <- readRDS(X_fn)


#Y <- read.csv(X_fn, check.names = FALSE)
#target_strings <- readLines(y_fn)
#stopifnot(nrow(Y) == length(target_strings))


# =========================================================
# helper: parse intervention targets
# =========================================================
parse_targets <- function(target_strings, colnames_Y) {
  targets <- lapply(target_strings, function(x) {
    x <- trimws(x)
    if (x == "control") {
      integer(0)
    } else {
      toks <- unlist(strsplit(x, "\\s+"))
      idx <- match(toks, colnames_Y)
      if (any(is.na(idx))) {
        stop(sprintf(
          "Target(s) not found in Y columns: %s",
          paste(toks[is.na(idx)], collapse = ", ")
        ))
      }
      idx
    }
  })
  
  target_labels <- sapply(targets, function(x) paste0(sort(x), collapse = ","))
  unique_target_keys <- unique(target_labels)
  
  unique_targets <- lapply(unique_target_keys, function(k) {
    if (k == "") integer(0) else as.integer(strsplit(k, ",")[[1]])
  })
  
  target.index <- match(target_labels, unique_target_keys)
  
  list(
    unique_targets = unique_targets,
    target.index = target.index
  )
}

# =========================================================
# helper: fit GIES on train split
# =========================================================
fit_gies_on_train <- function(Y_train, y_train, lam_base) {
  n_train <- nrow(Y_train)
  lam <- lam_base * log(n_train)
  
  parsed_train <- parse_targets(y_train, colnames(Y_train))
  
  score_obj_train <- new(
    "GaussL0penIntScore",
    data = Y_train,
    targets = parsed_train$unique_targets,
    target.index = parsed_train$target.index,
    lambda = lam
  )
  
  t0 <- Sys.time()
  gies_res <- gies(score_obj_train)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  list(
    gies_res = gies_res,
    lambda = lam,
    elapsed = elapsed
  )
}

# =========================================================
# helper: held-out score on test split
# IMPORTANT:
# fit graph on train, then evaluate that DAG on held-out test
# using a new score object built from the test data
# =========================================================
heldout_global_score <- function(Y_test, y_test, dag_repr, lam_base) {
  n_test <- nrow(Y_test)
  lam_test <- lam_base * log(n_test)
  
  parsed_test <- parse_targets(y_test, colnames(Y_test))
  
  score_obj_test <- new(
    "GaussL0penIntScore",
    data = Y_test,
    targets = parsed_test$unique_targets,
    target.index = parsed_test$target.index,
    lambda = lam_test
  )
  
  score_obj_test$global.score(dag_repr)
}

# =========================================================
# make separate tuning folds from FULL dataset
# seed 456 = tuning CV
# =========================================================

base_dir <- "/data/path/"

datasets <- list()

# 150d, 250d, 500d
for (d in c("150d", "250d", "500d")) {
  for (graph_type in c("random", "sf")) {
    out_dir <- file.path(base_dir, d, graph_type, "42")
    
    datasets[[length(datasets) + 1]] <- list(
      dim = d,
      setting = NA,
      graph_type = graph_type,
      out_dir = out_dir,
      X_fn = file.path(out_dir, "Y_matrix.csv"),
      y_fn = file.path(out_dir, "targets.txt")
    )
  }
}

# 50d special structure
for (setting in c("5", "15", "25", "50", "75", "100")) {
  for (graph_type in c("random", "sf")) {
    out_dir <- file.path(base_dir, "50d", setting, graph_type, "42")
    
    datasets[[length(datasets) + 1]] <- list(
      dim = "50d",
      setting = setting,
      graph_type = graph_type,
      out_dir = out_dir,
      X_fn = file.path(out_dir, "Y_matrix.csv"),
      y_fn = file.path(out_dir, "targets.txt")
    )
  }
}

length(datasets)

for (ds in datasets) {
  cat("\n====================================\n")
  cat("Running:", ds$dim, ds$setting, ds$graph_type, "\n")
  cat("Path:", ds$out_dir, "\n")
  cat("====================================\n")
  
  Y <- read.csv(ds$X_fn, check.names = FALSE)
  target_strings <- readLines(ds$y_fn)
  stopifnot(nrow(Y) == length(target_strings))
  
  idx_ctrl <- which(target_strings == "control")
  idx_pert <- which(target_strings != "control")
  
  set.seed(456)
  
  folds_pert <- caret::createFolds(
    factor(target_strings[idx_pert]),
    k = 5,
    returnTrain = FALSE
  )
  
  folds_ctrl <- caret::createFolds(
    seq_along(idx_ctrl),
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
    
    fold_scores <- numeric(5)
    fold_times <- numeric(5)
    
    for (fold in 1:5) {
      te_idx_pert <- idx_pert[folds_pert[[fold]]]
      te_idx_ctrl <- idx_ctrl[folds_ctrl[[fold]]]
      te_idx <- sort(c(te_idx_pert, te_idx_ctrl))
      
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
      
      fit <- fit_gies_on_train(Y_tr, y_tr, lam_base)
      
      score_val <- heldout_global_score(
        Y_test = Y_te,
        y_test = y_te,
        dag_repr = fit$gies_res$repr,
        lam_base = lam_base
      )
      
      fold_scores[fold] <- score_val
      fold_times[fold] <- fit$elapsed
      
      cat(sprintf(
        "   done | heldout_global_score=%.6f | gies_time=%.1fs\n",
        score_val, fit$elapsed
      ))
      flush.console()
    }
    
    mean_score <- mean(fold_scores)
    sd_score <- sd(fold_scores)
    
    cat(sprintf(
      ">>> lam_base=%.6f | mean heldout score=%.6f | sd=%.6f\n",
      lam_base, mean_score, sd_score
    ))
    flush.console()
    
    all_results[[g]] <- data.frame(
      lam_base = lam_base,
      fold1_score = fold_scores[1],
      fold2_score = fold_scores[2],
      fold3_score = fold_scores[3],
      fold4_score = fold_scores[4],
      fold5_score = fold_scores[5],
      mean_heldout_score = mean_score,
      sd_heldout_score = sd_score,
      mean_gies_time_sec = mean(fold_times)
    )
  }
  
  results_df <- do.call(rbind, all_results)
  results_df <- results_df[order(-results_df$mean_heldout_score), ]
  
  cat("\n================ FINAL TUNING RESULTS ================\n")
  print(results_df)
  
  best_lam_base <- results_df$lam_base[1]
  cat(sprintf("\nBest global lam_base = %.6f\n", best_lam_base))
}
