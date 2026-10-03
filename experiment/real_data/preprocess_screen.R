#!/usr/bin/env Rscript
#
# Preprocess one K562 Perturb-seq screen (essential or GWPS) for the real-data
# analysis (paper Section 5 / Appendix D.3). This replaces the original
# essential-only script, experiment/simulation/data_hdf5.R (removed; see git
# history), with the paths turned into arguments and the screen name in the
# output files. The processing is unchanged, and the output is byte-identical:
#
#   1. log1p(counts / core_scale_factor) on Replogle et al.'s raw counts;
#   2. Seurat S and G2M cell-cycle scores (100 expression-matched control genes
#      per cell-cycle gene from 24 mean-expression bins, seed 1234) regressed
#      out of each selected gene by least squares over all cells;
#   3. per GEM group, each gene z-scored against that group's non-targeting
#      cells;
#   4. only non-targeting cells and cells targeting a selected gene are kept.
#
# It then writes the IBCD inputs: the full data, and 5 cross-validation folds
# (caret::createFolds, stratified by target, seed 123), each as cells x genes
# with a final 'target' column ("control" for non-targeting cells).
#
# Usage:
#
#   Rscript preprocess_screen.R --raw_h5ad K562_essential_raw_singlecell_01.h5ad \
#       --genes common_genes.csv --screen essential --out_dir essential
#
# --genes is a CSV with a column common_genes whose entries contain Ensembl
# IDs (experiment/real_data/common_genes.csv, the published 521 genes, or the
# common_genes.csv written by select_genes.R).
#
# Writes, in out_dir:
#   X_ccc_norm_cellsxgenes_<screen>.rds, targets_clean_<screen>.txt
#   input/Y_matrix_<screen>_train.csv                      all cells
#   input/Y_matrix_<screen>_{train,test}_fold{1..5}.csv    CV folds

suppressMessages({
  library(hdf5r)
  library(tidyverse)
  library(matrixStats)
  library(caret)
})

get_arg <- function(args, key, default = NULL) {
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}
args <- commandArgs(trailingOnly = TRUE)
raw_h5ad <- get_arg(args, "--raw_h5ad")
gene_csv <- get_arg(args, "--genes")
screen   <- get_arg(args, "--screen")
out_dir  <- get_arg(args, "--out_dir", screen)
if (is.null(raw_h5ad) || is.null(gene_csv) || is.null(screen)) stop("need --raw_h5ad, --genes and --screen")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ts <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

ensg_ids <- readr::read_csv(gene_csv, show_col_types = FALSE)$common_genes |>
  str_extract("ENSG[0-9]{11}")
ts(screen, ": ", length(ensg_ids), " genes from ", gene_csv)

h <- H5File$new(raw_h5ad, "r")

cats_obs <- lapply(names(h[["obs/__categories"]]), function(nm) {
  vals <- h[["obs/__categories"]][[nm]][]
  set_names(as.character(0:(length(vals)-1)), vals)
}) |> set_names(names(h[["obs/__categories"]]))

# obs & var
obs <- tibble(
  cell_barcode    = h[["obs/cell_barcode"]][],
  core_scale      = h[["obs/core_scale_factor"]][],
  gem_group       = h[["obs/gem_group"]][],
  gene            = h[["obs/gene"]][],
  gene_id         = h[["obs/gene_id"]][],
  gene_transcript = h[["obs/gene_transcript"]][]
) |>
  mutate(
    gene            = fct_recode(factor(gene),            !!!cats_obs$gene)            |> as.character(),
    gene_id         = fct_recode(factor(gene_id),         !!!cats_obs$gene_id)         |> as.character(),
    gene_transcript = fct_recode(factor(gene_transcript), !!!cats_obs$gene_transcript) |> as.character()
  )

cats_var <- lapply(names(h[["var/__categories"]]), function(nm) {
  vals <- h[["var/__categories"]][[nm]][]
  set_names(as.character(0:(length(vals)-1)), vals)
}) |> set_names(names(h[["var/__categories"]]))

var <- tibble(
  gene_id   = h[["var/gene_id"]][],
  gene_name = h[["var/gene_name"]][]
) |>
  mutate(gene_name = fct_recode(factor(gene_name), !!!cats_var$gene_name) |> as.character())

G <- nrow(var); C <- nrow(obs)
ts("dims genes×cells: ", G, "×", C)
stopifnot(all(dim(h[["X"]]) == c(G, C)))

# ---- helpers ----
# read a block (arbitrary row & contiguous col range), log1p(core-scale)
read_rows_cols <- function(row_idx, c0, c1) {
  X <- h[["X"]][row_idx, c0:c1, drop = FALSE]
  cs <- obs$core_scale[c0:c1]
  log1p(sweep(X, 2, cs, "/"))
}

# ---- cell cycle gene sets (Seurat S/G2M) mapped to ENSG ----
cc_genes <- list(
  s   = c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL",
          "PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP",
          "CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN",
          "DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"),
  g2m = c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B",
          "MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1",
          "KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C",
          "KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA",
          "PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
)
sym2ensg <- distinct(var, gene_name, gene_id)
S_idx   <- match(sym2ensg$gene_id[sym2ensg$gene_name %in% cc_genes$s],   var$gene_id) |> (\(x) x[!is.na(x)])()
G2M_idx <- match(sym2ensg$gene_id[sym2ensg$gene_name %in% cc_genes$g2m], var$gene_id) |> (\(x) x[!is.na(x)])()
ts("cell-cycle genes found: S ", length(S_idx), ", G2M ", length(G2M_idx))

# ---- (A) per-gene means over all cells (for binning) ----
ts("per-gene means (binning)…")
gene_sums <- numeric(G)
c_bs_means <- 5000L    # safe block for passes that touch many genes
for (c0 in seq(1L, C, by = c_bs_means)) {
  c1 <- min(c0 + c_bs_means - 1L, C)
  Yc <- h[["X"]][, c0:c1, drop = FALSE]
  Yc <- log1p(sweep(Yc, 2, obs$core_scale[c0:c1], "/"))
  gene_sums <- gene_sums + rowSums(Yc)
  rm(Yc); gc()
  if (((c0 - 1L) / c_bs_means) %% 20 == 0) ts("  cells ", c1, "/", C)
}
avg <- gene_sums / C

# ---- (B) pick matched control genes per CC set ----
bins <- cut_number(avg, 24); names(bins) <- var$gene_id
pick_ctrl <- function(gene_idx_vec, n_ctrl = 100) {
  gids <- var$gene_id[gene_idx_vec]
  ctrl <- unique(unlist(lapply(gids, function(g) {
    b <- bins[g]; if (is.na(b)) character(0) else sample(names(bins)[bins == b], n_ctrl)
  })))
  match(ctrl, var$gene_id) |> (\(x) x[!is.na(x)])()
}
set.seed(1234)
S_ctrl_idx   <- pick_ctrl(S_idx)
G2M_ctrl_idx <- pick_ctrl(G2M_idx)

# ---- (C) CC module scores for every cell (read only needed rows) ----
ts("CC module scores…")
need_rows <- sort(unique(c(S_idx, S_ctrl_idx, G2M_idx, G2M_ctrl_idx)))
S_score   <- numeric(C)
G2M_score <- numeric(C)
for (c0 in seq(1L, C, by = c_bs_means)) {
  c1 <- min(c0 + c_bs_means - 1L, C)
  Yc_all <- read_rows_cols(need_rows, c0, c1)
  # map back
  map_rows <- function(idx) match(idx, need_rows)
  S_mat    <- Yc_all[map_rows(S_idx), , drop = FALSE]
  S_ctrl   <- Yc_all[map_rows(S_ctrl_idx), , drop = FALSE]
  G2M_mat  <- Yc_all[map_rows(G2M_idx), , drop = FALSE]
  G2M_ctrl <- Yc_all[map_rows(G2M_ctrl_idx), , drop = FALSE]
  S_score[c0:c1]   <- colMeans(S_mat)   - colMeans(S_ctrl)
  G2M_score[c0:c1] <- colMeans(G2M_mat) - colMeans(G2M_ctrl)
  rm(Yc_all, S_mat, S_ctrl, G2M_mat, G2M_ctrl); gc()
}

# ---- (D) CC design & (X'X)^{-1} ----
CC <- cbind(S_score = S_score, G2M_score = G2M_score, Intercept = 1)
XtX_inv <- solve(crossprod(CC))   # 3×3

# ---- (E) which genes/cells we keep ----
keep_genes   <- intersect(ensg_ids, var$gene_id)
keep_gene_ix <- match(keep_genes, var$gene_id) |> (\(x) x[!is.na(x)])()
cells_filter <- obs$cell_barcode[obs$gene_id == "non-targeting" | obs$gene_id %in% ensg_ids]
keep_cell_ix <- match(cells_filter, obs$cell_barcode) |> (\(x) x[!is.na(x)])()
ts("keeping ", length(keep_gene_ix), " of ", length(ensg_ids), " genes and ",
   length(keep_cell_ix), " cells (", sum(obs$gene_id == "non-targeting"), " non-targeting)")
stopifnot(!anyDuplicated(obs$cell_barcode[keep_cell_ix]))

# pre-allocate final X (cells × genes)
X <- matrix(NA_real_, nrow = length(keep_cell_ix), ncol = length(keep_gene_ix),
            dimnames = list(obs$cell_barcode[keep_cell_ix], var$gene_id[keep_gene_ix]))
row_map <- setNames(seq_along(keep_cell_ix), obs$cell_barcode[keep_cell_ix])

# ---- (F) estimate CC betas for the kept genes (streaming over cells) ----
ts("estimating betas for ", length(keep_gene_ix), " genes…")
g_bs <- 300L       # gene batch
c_bs <- 20000L     # cell batch (adjust if RAM is tight)
betas <- vector("list", length(keep_gene_ix))  # each a length-3 numeric

for (gi0 in seq(1L, length(keep_gene_ix), by = g_bs)) {
  gi1  <- min(gi0 + g_bs - 1L, length(keep_gene_ix))
  rows <- keep_gene_ix[gi0:gi1]
  glen <- length(rows)

  S <- matrix(0, nrow = ncol(CC), ncol = glen)  # 3×glen
  for (c0 in seq(1L, C, by = c_bs)) {
    c1  <- min(c0 + c_bs - 1L, C)
    Ygc <- read_rows_cols(rows, c0, c1)         # glen×cells_blk
    CCc <- CC[c0:c1, , drop = FALSE]            # cells_blk×3
    S   <- S + crossprod(CCc, t(Ygc))           # 3×glen
    rm(Ygc, CCc); gc()
  }
  B <- XtX_inv %*% S  # 3×glen
  for (k in seq_len(glen)) betas[[gi0 + k - 1L]] <- as.numeric(B[, k])
  rm(S, B); gc()
}

# ---- (G) GEM-normalize per gem_group and fill X directly ----
.pct <- function(done, total) sprintf("%5.1f%%", 100 * done / max(total, 1))
.eta <- function(t0, done, total) {
  if (done <= 0) return("ETA --:--")
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rem <- elapsed * (total / done - 1)
  sprintf("ETA %02d:%02d", floor(rem / 60), round(rem %% 60))
}

ts("GEM-normalize into X…")

Gkeep <- length(keep_gene_ix)
gem_levels <- sort(unique(obs$gem_group))
g_total <- length(gem_levels)
g_start_time <- Sys.time()

for (gi in seq_along(gem_levels)) {
  g <- gem_levels[gi]
  all_cols <- which(obs$gem_group == g)
  ntc_cols <- which(obs$gem_group == g & obs$gene_id == "non-targeting")
  keep_in_g <- which(obs$gem_group == g & obs$cell_barcode %in% names(row_map))

  ts(sprintf("GEM %s  (%d/%d): all=%d, NTC=%d, kept=%d",
             g, gi, g_total, length(all_cols), length(ntc_cols), length(keep_in_g)))

  if (length(all_cols) == 0L || length(ntc_cols) < 2L) {
    ts(sprintf("GEM %s: skipped (no cells or <2 NTC).", g))
    next
  }

  # ---- pass 1: NTC mean/sd per gene ----
  mu  <- numeric(Gkeep)
  vv  <- numeric(Gkeep)
  ntc_total <- length(ntc_cols)
  pass1_t0  <- Sys.time()

  for (gi0 in seq(1L, Gkeep, by = g_bs)) {
    gi1  <- min(gi0 + g_bs - 1L, Gkeep)
    rows <- keep_gene_ix[gi0:gi1]
    glen <- length(rows)
    B    <- do.call(cbind, betas[gi0:gi1])     # 3×glen

    sumR  <- numeric(glen)
    sumR2 <- numeric(glen)

    for (c0 in seq(1L, ntc_total, by = c_bs)) {
      idx <- ntc_cols[c0:min(c0 + c_bs - 1L, ntc_total)]
      Ygc <- read_rows_cols(rows, min(idx), max(idx))[
        , (idx - min(idx) + 1L), drop = FALSE]     # glen×cells_blk
      CCc <- CC[idx, , drop = FALSE]                      # cells_blk×3
      R   <- Ygc - t(CCc %*% B)

      sumR  <- sumR  + rowSums(R)
      sumR2 <- sumR2 + rowSums(R * R)
      rm(Ygc, CCc, R); gc(FALSE)
    }

    mu_batch      <- sumR / ntc_total
    var_batch     <- pmax(sumR2 / ntc_total - mu_batch^2, 1e-12)
    mu[gi0:gi1]   <- mu_batch
    vv[gi0:gi1]   <- var_batch
    rm(B, sumR, sumR2, mu_batch, var_batch); gc(FALSE)
  }
  sdv <- sqrt(vv)

  # ---- pass 2: normalize kept cells and write to X ----
  k_total <- length(keep_in_g)
  if (k_total == 0L) {
    ts(sprintf("GEM %s: no kept cells; skipping write.", g))
    next
  }

  for (c0 in seq(1L, k_total, by = c_bs)) {
    idx      <- keep_in_g[c0:min(c0 + c_bs - 1L, k_total)]
    CCc      <- CC[idx, , drop = FALSE]
    out_rows <- unname(row_map[obs$cell_barcode[idx]])

    for (gi0 in seq(1L, Gkeep, by = g_bs)) {
      gi1  <- min(gi0 + g_bs - 1L, Gkeep)
      rows <- keep_gene_ix[gi0:gi1]
      B    <- do.call(cbind, betas[gi0:gi1])               # 3×glen

      Ygc <- read_rows_cols(rows, min(idx), max(idx))[
        , (idx - min(idx) + 1L), drop = FALSE]      # glen×cells_blk
      R   <- Ygc - t(CCc %*% B)
      Z   <- (R - mu[gi0:gi1]) / sdv[gi0:gi1]
      X[out_rows, gi0:gi1] <- t(Z)
      rm(Ygc, R, Z, B); gc(FALSE)
    }
    rm(CCc); gc(FALSE)
  }

  ts(sprintf("GEM %s: done. (%.1f sec)  overall %s  %s", g,
             as.numeric(difftime(Sys.time(), pass1_t0, units = "secs")),
             .pct(gi, g_total), .eta(g_start_time, gi, g_total)))
}
h$close_all()

# ---- targets & save ----
obs_keep <- obs[keep_cell_ix, ]
targets_clean <- ifelse(obs_keep$gene_id == "non-targeting", "control", obs_keep$gene_id)
stopifnot(nrow(X) == length(targets_clean))

# cells in a GEM group that was skipped (fewer than 2 NTC) were never filled
unfilled <- rowSums(is.na(X)) > 0
if (any(unfilled)) {
  ts("dropping ", sum(unfilled), " cells left unnormalised (GEM groups with < 2 NTC)")
  X <- X[!unfilled, , drop = FALSE]
  targets_clean <- targets_clean[!unfilled]
}

saveRDS(X, file.path(out_dir, sprintf("X_ccc_norm_cellsxgenes_%s.rds", screen)))
writeLines(targets_clean, file.path(out_dir, sprintf("targets_clean_%s.txt", screen)))
ts("dim(X) = ", paste(dim(X), collapse = "×"),
   " | targets: ", sum(targets_clean == "control"), " control, ",
   length(unique(targets_clean[targets_clean != "control"])), " perturbed genes")

# ---- IBCD inputs: full data and 5 CV folds ----
input_dir <- file.path(out_dir, "input")
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

write_input <- function(rows, name) {
  Y <- as.data.frame(X[rows, , drop = FALSE])
  Y$target <- targets_clean[rows]
  out_csv <- file.path(input_dir, sprintf("Y_matrix_%s_%s.csv", screen, name))
  write.csv(Y, out_csv, row.names = FALSE)
  ts(sprintf("%s: %d samples (ctrl=%d, pert=%d)", basename(out_csv), length(rows),
             sum(Y$target == "control"), sum(Y$target != "control")))
}

write_input(seq_len(nrow(X)), "train")

idx_ctrl <- which(targets_clean == "control")
idx_pert <- which(targets_clean != "control")

# CV performance stability
set.seed(123)
folds_pert <- caret::createFolds(factor(targets_clean[idx_pert]), k = 5, returnTrain = FALSE)
folds_ctrl <- caret::createFolds(seq_along(idx_ctrl), k = 5, returnTrain = FALSE)

for (i in seq_len(5)) {
  te_idx <- sort(c(idx_pert[folds_pert[[i]]], idx_ctrl[folds_ctrl[[i]]]))
  tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
  stopifnot(length(tr_idx) + length(te_idx) == nrow(X))
  write_input(tr_idx, sprintf("train_fold%d", i))
  write_input(te_idx, sprintf("test_fold%d", i))
}
ts("done")
