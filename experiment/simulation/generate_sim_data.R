#!/usr/bin/env Rscript
#
# Generate interventional simulation datasets for the IBCD experiments.
#
# Each (D, n_int, graph, seed) combination gets its own directory:
#
#     <out_dir>/<D>d/<n_int>/<graph>/<seed>/
#
# containing:
#
#     Y_matrix.csv   data matrix, no target column
#     targets.txt    one label per row, "control" or "V<j>"
#     G_matrix.csv   true direct effects
#     R_matrix.csv   true total effects
#     xi_true.csv    normalised squared control covariance
#     R_true.csv, SE_hat.csv, U_true.csv, V_true.csv   2SLS summaries
#
# and additionally:
#
#     Y_with_targets.csv   Y_matrix.csv with the targets joined on as a column
#
# ibcd.py needs the target column; the other methods take Y_matrix.csv and
# targets.txt separately, so both layouts are written rather than changing the
# existing one.
#
# Usage:
#
#   Rscript generate_sim_data.R --out_dir DIR [--dims 50,150,250,500]
#       [--graphs er,sf] [--n_int 100] [--seeds 42:51] [--shuffle false]
#       [--write_iv true]
#
#
# Node order
# ----------
# inspre's DAG generators emit a graph that is upper triangular in the order
# the variables are written, so column order is a topological order of the true
# graph. That is an artefact of the simulator and would never hold in real
# data. --shuffle true permutes the variables after generation so no method can
# benefit from the ordering, but it changes the datasets, which invalidates
# results computed on the unshuffled ones. Off by default for that reason.

suppressMessages({
  library(inspre)
  library(mc2d)
})

# ---------------------------------------------------------------- args ----

get_arg <- function(args, key, default = NULL) {
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}

parse_list <- function(s) trimws(strsplit(s, ",", fixed = TRUE)[[1]])

parse_seeds <- function(s) {
  s <- trimws(s)
  if (grepl(":", s, fixed = TRUE)) {
    p <- as.integer(strsplit(s, ":", fixed = TRUE)[[1]])
    seq(p[1], p[2])
  } else {
    as.integer(parse_list(s))
  }
}

parse_bool <- function(s) tolower(trimws(s)) %in% c("true", "t", "yes", "1")

args <- commandArgs(trailingOnly = TRUE)
out_dir <- get_arg(args, "--out_dir")
if (is.null(out_dir)) stop("--out_dir is required")

dims     <- as.integer(parse_list(get_arg(args, "--dims", "50,150,250,500")))
graphs   <- parse_list(get_arg(args, "--graphs", "er,sf"))
n_ints   <- as.integer(parse_list(get_arg(args, "--n_int", "100")))
seeds    <- parse_seeds(get_arg(args, "--seeds", "42:51"))
shuffle  <- parse_bool(get_arg(args, "--shuffle", "false"))
write_iv <- parse_bool(get_arg(args, "--write_iv", "true"))
int_beta <- as.numeric(get_arg(args, "--int_beta", "-2"))
v_edge   <- as.numeric(get_arg(args, "--v", "0.25"))

stopifnot(all(graphs %in% c("er", "sf")))

# Appendix D.1: p chosen per D so average degree is about 5 for each family.
P_GRID <- list(
  er = c("50" = 0.10,  "150" = 0.033, "250" = 0.020, "500" = 0.010),
  sf = c("50" = 0.066, "150" = 0.108, "250" = 0.124, "500" = 0.139)
)
# inspre's generate_network takes 'random' for Erdos-Renyi.
INSPRE_GRAPH <- c(er = "random", sf = "scalefree")

# ------------------------------------------------------------ helpers ----

multiple_iv_reg_UV <- function(target, .X, .targets) {
  inst_obs <- .targets == target
  control_obs <- .targets == "control"
  n_target <- sum(inst_obs)
  X_target <- .X[inst_obs, , drop = FALSE]
  X_control <- .X[control_obs, , drop = FALSE]
  this_feature <- which(colnames(.X) == target)
  X_exp_target <- X_target[, this_feature]
  X_exp_control <- X_control[, this_feature]
  beta_inst_obs <- sum(X_exp_target) / n_target
  beta_hat <- colSums(X_target) / colSums(X_target)[this_feature]
  resid <- rbind(X_target - outer(X_exp_target, beta_hat),
                 X_control - outer(X_exp_control, beta_hat))
  V_hat <- crossprod(resid) / nrow(resid)
  se_hat <- sqrt(diag(V_hat) / (n_target * beta_inst_obs^2))
  list(beta_hat = beta_hat, se_hat = se_hat,
       U_i = 1 / (n_target * beta_inst_obs^2), V_hat = V_hat)
}

#' Relabel the variables by a permutation, keeping Y, G, R and targets aligned.
#'
#' New variable k is old variable perm[k], so a sample intervened on old V_j
#' becomes V_{inv[j]} with inv = order(perm).
permute_dataset <- function(Y, G, R, targets, perm) {
  D <- length(perm)
  inv <- order(perm)
  Y <- Y[, perm, drop = FALSE]
  colnames(Y) <- paste0("V", seq_len(D))
  is_ctrl <- targets == "control"
  old_idx <- as.integer(sub("^V", "", targets[!is_ctrl]))
  targets[!is_ctrl] <- paste0("V", inv[old_idx])
  list(Y = Y, G = G[perm, perm, drop = FALSE], R = R[perm, perm, drop = FALSE],
       targets = targets)
}

run_one <- function(D, n_int, graph, seed, dir_out) {
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  p <- P_GRID[[graph]][[as.character(D)]]
  if (is.null(p) || is.na(p)) {
    stop(sprintf("no p defined for graph=%s D=%d; add it to P_GRID", graph, D))
  }
  message(sprintf("[D=%d | n_int=%d | %s (p=%.3f) | seed=%d] -> %s",
                  D, n_int, graph, p, seed, dir_out))
  set.seed(seed)

  ds <- generate_dataset(
    D = D, N_cont = D * n_int, N_int = n_int, int_beta = int_beta,
    graph = INSPRE_GRAPH[[graph]],     # passed explicitly; never defaulted
    v = v_edge, p = p, DAG = TRUE, C = 0
  )
  Y <- ds$Y; G <- ds$G; R <- ds$R; targets <- ds$targets

  if (shuffle) {
    pd <- permute_dataset(Y, G, R, targets, sample.int(D))
    Y <- pd$Y; G <- pd$G; R <- pd$R; targets <- pd$targets
  }

  Y_C <- Y[targets == "control", , drop = FALSE]
  xi <- (crossprod(Y_C) / nrow(Y_C))^2
  diag(xi) <- 0
  xi <- xi / sqrt(sum(xi^2))

  writeLines(targets, file.path(dir_out, "targets.txt"))
  write.csv(Y,  file.path(dir_out, "Y_matrix.csv"), row.names = FALSE)
  write.csv(G,  file.path(dir_out, "G_matrix.csv"), row.names = FALSE)
  write.csv(R,  file.path(dir_out, "R_matrix.csv"), row.names = FALSE)
  write.csv(xi, file.path(dir_out, "xi_true.csv"),  row.names = FALSE)

  # ibcd.py reads a single CSV carrying the targets as a column
  Y_ibcd <- as.data.frame(Y)
  Y_ibcd$target <- targets
  write.csv(Y_ibcd, file.path(dir_out, "Y_with_targets.csv"), row.names = FALSE)

  if (write_iv) {
    genes <- paste0("V", seq_len(D))
    U_diag <- numeric(D); V_sum <- matrix(0, D, D)
    R_hat <- matrix(0, D, D); SE_hat <- matrix(0, D, D)
    for (i in seq_len(D)) {
      res <- multiple_iv_reg_UV(genes[i], Y, targets)
      U_diag[i] <- res$U_i
      V_sum <- V_sum + res$V_hat
      R_hat[i, ] <- res$beta_hat
      SE_hat[i, ] <- res$se_hat
    }
    write.csv(R_hat,        file.path(dir_out, "R_true.csv"),  row.names = FALSE)
    write.csv(SE_hat,       file.path(dir_out, "SE_hat.csv"),  row.names = FALSE)
    write.csv(diag(U_diag), file.path(dir_out, "U_true.csv"),  row.names = FALSE)
    write.csv(V_sum / D,    file.path(dir_out, "V_true.csv"),  row.names = FALSE)
  }
}

# --------------------------------------------------------------- main ----

message(sprintf("out_dir=%s | dims=%s | graphs=%s | n_int=%s | seeds=%s | shuffle=%s",
                out_dir, paste(dims, collapse = ","), paste(graphs, collapse = ","),
                paste(n_ints, collapse = ","),
                paste(range(seeds), collapse = ":"), shuffle))

for (D in dims) {
  for (n_int in n_ints) {
    for (graph in graphs) {
      for (seed in seeds) {
        run_one(D, n_int, graph, seed,
                file.path(out_dir, paste0(D, "d"), as.character(n_int),
                          graph, as.character(seed)))
      }
    }
  }
}

message("Done.")
