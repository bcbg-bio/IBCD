#!/usr/bin/env Rscript
#
# SID for one estimated graph, computed exactly as experiment/simulation/
# metrics.R does: both matrices thresholded at |G| > eps and passed straight
# to SID::structIntervDist, with no resolution of cycles. Split out so a
# cluster array can run one graph per task, since a single SID call at
# D = 500 can take hours when the estimate has many cycles.
#
# Usage:
#
#   Rscript sid_one.R --g_true G_matrix.csv --g_est G.csv --out sid.csv
#       [--eps 0.075] [--label "D=500,n_int=100,graph=sf,seed=42"]
#
# Writes one CSV row: label, sid, seconds. The estimate may carry a row-name
# column (as ibcd.py's G.csv does); it is detected and dropped.

suppressMessages(library(SID))

get_arg <- function(args, key, default = NULL) {
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}
args <- commandArgs(trailingOnly = TRUE)
g_true <- get_arg(args, "--g_true")
g_est  <- get_arg(args, "--g_est")
out    <- get_arg(args, "--out")
eps    <- as.numeric(get_arg(args, "--eps", "0.075"))
label  <- get_arg(args, "--label", "")
if (is.null(g_true) || is.null(g_est) || is.null(out)) stop("need --g_true, --g_est and --out")

read_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  if (!is.numeric(df[[1]]) || ncol(df) == nrow(df) + 1) df <- read.csv(path, row.names = 1, check.names = FALSE)
  as.matrix(df)
}
G_true_mat <- read_matrix(g_true)
G_est_mat  <- read_matrix(g_est)
stopifnot(all(dim(G_true_mat) == dim(G_est_mat)))

A_true <- 1 * (abs(G_true_mat) > eps)
A_est  <- 1 * (abs(G_est_mat)  > eps)

t0 <- Sys.time()
sid_val <- structIntervDist(A_true, A_est)$sid
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

write.csv(data.frame(label = label, sid = sid_val, seconds = secs, status = "ok"),
          out, row.names = FALSE)
cat(sprintf("%s SID %s (%.1f s)\n", label, sid_val, secs))
