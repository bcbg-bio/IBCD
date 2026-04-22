library(hdf5r)
library(inspre)
library(tidyverse)
library(forcats)
library(stringr)
library(Matrix)
library(matrixStats) 
library(scales)
library(caret)

raw_h5ad <- "K562_essential_raw_singlecell_01.h5ad"
gene_csv <- "test.csv"
out_dir  <- "output"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ts <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

ensg_ids <- readr::read_csv(gene_csv, show_col_types = FALSE)$common_genes |>
  str_extract("ENSG[0-9]{11}")

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

# pre-allocate final X (cells × genes)
X <- matrix(NA_real_, nrow = length(keep_cell_ix), ncol = length(keep_gene_ix),
            dimnames = list(obs$cell_barcode[keep_cell_ix], var$gene_id[keep_gene_ix]))
row_map <- setNames(seq_along(keep_cell_ix), obs$cell_barcode[keep_cell_ix])

# ---- (F) estimate betas for the 521 genes (streaming over cells) ----
ts("estimating betas for 521 genes…")
g_bs <- 300L       # gene batch
c_bs <- 20000L     # cell batch for 521 genes (adjust if RAM is tight)
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
# ---------- VERBOSE GEM-NORMALIZE (drop-in replacement) ----------
ts <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
.pct <- function(done, total) sprintf("%5.1f%%", 100 * done / max(total, 1))
.eta <- function(t0, done, total) {
  if (done <= 0) return("ETA --:--")
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rem <- elapsed * (total / done - 1)
  sprintf("ETA %02d:%02d", floor(rem / 60), round(rem %% 60))
}

ts("GEM-normalize into X…")

g_bs <- 300L        # gene batch (same as above)
c_bs <- 20000L      # cell batch (same as above)
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
  ts(sprintf("GEM %s: NTC mean/sd (genes=%d, ntc_cells=%d)…",
             g, Gkeep, length(ntc_cols)))
  mu  <- numeric(Gkeep)
  vv  <- numeric(Gkeep)
  
  ntc_total <- length(ntc_cols)
  ntc_done  <- 0L
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
      
      ntc_done <- ntc_done + ncol(R)
      ts(sprintf("  NTC %7d/%7d  (%s)  genes %4d-%-4d/%-4d   %s",
                 ntc_done, ntc_total, .pct(ntc_done, ntc_total),
                 gi0, gi1, Gkeep, .eta(pass1_t0, ntc_done, ntc_total)))
      
      rm(Ygc, CCc, R); gc(FALSE)
    }
    
    mu_batch      <- sumR / ntc_total
    var_batch     <- pmax(sumR2 / ntc_total - mu_batch^2, 1e-12)
    mu[gi0:gi1]   <- mu_batch
    vv[gi0:gi1]   <- var_batch
    
    ts(sprintf("  genes %4d-%-4d/%-4d NTC stats committed.",
               gi0, gi1, Gkeep))
    rm(B, sumR, sumR2, mu_batch, var_batch); gc(FALSE)
  }
  sdv <- sqrt(vv)
  ts(sprintf("GEM %s: NTC mean/sd done. (%.1f sec)",
             g, as.numeric(difftime(Sys.time(), pass1_t0, units = "secs"))))
  
  # ---- pass 2: normalize kept cells and write to X ----
  k_total <- length(keep_in_g)
  if (k_total == 0L) {
    ts(sprintf("GEM %s: no kept cells; skipping write.", g))
    next
  }
  ts(sprintf("GEM %s: normalize kept cells (%d)…", g, k_total))
  
  pass2_t0 <- Sys.time()
  k_done   <- 0L
  for (c0 in seq(1L, k_total, by = c_bs)) {
    idx      <- keep_in_g[c0:min(c0 + c_bs - 1L, k_total)]
    CCc      <- CC[idx, , drop = FALSE]
    out_rows <- unname(row_map[obs$cell_barcode[idx]])
    
    for (gi0 in seq(1L, Gkeep, by = g_bs)) {
      gi1  <- min(gi0 + g_bs - 1L, Gkeep)
      rows <- keep_gene_ix[gi0:gi1]
      glen <- length(rows)
      B    <- do.call(cbind, betas[gi0:gi1])               # 3×glen
      
      Ygc <- read_rows_cols(rows, min(idx), max(idx))[
        , (idx - min(idx) + 1L), drop = FALSE]      # glen×cells_blk
      R   <- Ygc - t(CCc %*% B)
      Z   <- (R - mu[gi0:gi1]) / sdv[gi0:gi1]
      X[out_rows, gi0:gi1] <- t(Z)
      
      ts(sprintf("  cells %7d/%7d  (%s)  genes %4d-%-4d/%-4d",
                 min(c0 + c_bs - 1L, k_total), k_total, .pct(min(c0 + c_bs - 1L, k_total), k_total),
                 gi0, gi1, Gkeep))
      
      rm(Ygc, R, Z, B); gc(FALSE)
    }
    
    k_done <- min(c0 + c_bs - 1L, k_total)
    ts(sprintf("  wrote cells %7d/%7d  (%s)  %s",
               k_done, k_total, .pct(k_done, k_total),
               .eta(pass2_t0, k_done, k_total)))
    rm(CCc); gc(FALSE)
  }
  
  ts(sprintf("GEM %s: done. (%.1f sec)", g,
             as.numeric(difftime(Sys.time(), pass2_t0, units = "secs"))))
  ts(sprintf("Overall progress %d/%d gems  (%s)  %s",
             gi, g_total, .pct(gi, g_total), .eta(g_start_time, gi, g_total)))
}
# ---------- END VERBOSE GEM-NORMALIZE ----------
h$close_all()

# ---- targets & save ----
obs_keep <- obs[keep_cell_ix, ]
targets_clean <- ifelse(obs_keep$gene_id == "non-targeting", "control", obs_keep$gene_id)

saveRDS(X, file.path(out_dir, "X_521_ccc_norm_cellsxgenes_essential.rds")) 
writeLines(targets_clean, file.path(out_dir, "targets_clean_essential.txt"))

ts("Done. dim(X) = ", paste(dim(X), collapse = "×"),
   " | length(targets_clean) = ", length(targets_clean))
stopifnot(nrow(X) == length(targets_clean))


dim(X)        # cells × ~521
length(targets_clean)
stopifnot(nrow(X) == length(targets_clean))

################################
################################

X_fn <- file.path(out_dir, "X_521_ccc_norm_cellsxgenes_essential.rds")
y_fn <- file.path(out_dir, "targets_clean_essential.txt")

# ---- load ----
X <- readRDS(X_fn)                    # cells × genes
targets_clean <- readLines(y_fn)      # same length as nrow(X)

stopifnot(nrow(X) == length(targets_clean))

# ---- split indices ----
idx_ctrl <- which(targets_clean == "control")
idx_pert <- which(targets_clean != "control")

# CV performance stability
set.seed(123)

# CV hyperparmeter tune
#set.seed(456)

folds_pert <- caret::createFolds(
  factor(targets_clean[idx_pert]),
  k = 5,
  returnTrain = FALSE
)

folds_ctrl <- caret::createFolds(
  seq_along(idx_ctrl),
  k = 5,
  returnTrain = FALSE
)

fold_dir <- "cv_data"
dir.create(fold_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(5)) {
  
  te_idx_pert <- idx_pert[folds_pert[[i]]]
  te_idx_ctrl <- idx_ctrl[folds_ctrl[[i]]]
  te_idx <- sort(c(te_idx_pert, te_idx_ctrl))
  
  tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
  
  X_tr <- X[tr_idx, , drop = FALSE]
  y_tr <- targets_clean[tr_idx]
  
  X_te <- X[te_idx, , drop = FALSE]
  y_te <- targets_clean[te_idx]
  
  stopifnot(nrow(X_tr) == length(y_tr))
  stopifnot(nrow(X_te) == length(y_te))
  stopifnot(length(intersect(tr_idx, te_idx)) == 0)
  stopifnot(nrow(X_tr) + nrow(X_te) == nrow(X))
  
  y_tr_out <- file.path(fold_dir, sprintf("Y_matrix_essential_train_fold%d.csv", i))
  t_tr_out <- file.path(fold_dir, sprintf("targets_essential_train_fold%d.txt", i))
  
  y_te_out <- file.path(fold_dir, sprintf("Y_matrix_essential_test_fold%d.csv", i))
  t_te_out <- file.path(fold_dir, sprintf("targets_essential_test_fold%d.txt", i))
  
  write.csv(X_tr, y_tr_out, row.names = FALSE)
  writeLines(y_tr, t_tr_out)
  
  write.csv(X_te, y_te_out, row.names = FALSE)
  writeLines(y_te, t_te_out)
  
  cat(sprintf(
    paste0(
      "   TRAIN: %d samples (ctrl=%d, pert=%d)\n",
      "   TEST : %d samples (ctrl=%d, pert=%d)\n"
    ),
    i,
    length(y_tr), sum(y_tr == "control"), sum(y_tr != "control"),
    length(y_te), sum(y_te == "control"), sum(y_te != "control")
  ))
}

out_dir  <- file.path(fold_dir, "input")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (i in seq_len(5)) {
  
  for (split in c("train", "test")) {
    
    y_csv <- file.path(fold_dir, sprintf("Y_matrix_essential_%s_fold%d.csv", split, i))
    t_txt <- file.path(fold_dir, sprintf("targets_essential_%s_fold%d.txt", split, i))

    out_csv <- file.path(out_dir, sprintf("Y_matrix_essential_%s_fold%d.csv", split, i))
    
    Y <- read.csv(y_csv, check.names = FALSE)
    targets <- readLines(t_txt)
    
    stopifnot(nrow(Y) == length(targets))
    Y$target <- targets
    
    write.csv(Y, out_csv, row.names = FALSE)
    
    cat(sprintf(
      "Created: %s from %s + %s\n",
      basename(out_csv), basename(y_csv), basename(t_txt)
    ))
  }
}


