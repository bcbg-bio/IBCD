#!/usr/bin/env Rscript
#
# Select the genes for the K562 Perturb-seq analysis (paper Appendix D.3).
#
# A gene is kept if, in BOTH the essential and the genome-wide (GWPS) screens,
#   - its best guide lowers the gene's own normalised expression by more than
#     0.75 standard deviations (inst_beta < -0.75), and
#   - more than 50 cells carry that guide (n > n_ntc + 50, where n counts the
#     non-targeting control cells plus the guide's cells),
# measured on Replogle et al.'s normalised single-cell data, i.e. before the
# cell-cycle regression. The paper reports 521 genes.
#
# The guide-effect calculation is the one in inspre's run_gwps_analysis.R:
# inspre::calc_inst_effect_h5X per guide, the guide with the smallest
# FDR-adjusted one-sided p-value kept per target. calc_inst_effect_h5X and
# parse_hdf5_df are copied here unchanged from inspre 1.0.1 so the cluster
# needs only hdf5r and dplyr, not inspre.
#
# Usage:
#
#   Rscript select_genes.R --essential K562_essential_normalized_singlecell_01.h5ad \
#       --gwps K562_gwps_normalized_singlecell_01.h5ad --out_dir genes \
#       [--max_beta -0.75] [--min_cells 50]
#
# Writes, in out_dir:
#   guide_effects_essential.csv, guide_effects_gwps.csv   best guide per target
#   kept_essential.csv, kept_gwps.csv                     targets passing both filters
#   common_genes.csv                                      column common_genes: the
#                                                         intersection, the input
#                                                         preprocess_screen.R reads

suppressMessages({library(hdf5r); library(dplyr)})

get_arg <- function(args, key, default = NULL) {
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}
args <- commandArgs(trailingOnly = TRUE)
files <- c(essential = get_arg(args, "--essential"), gwps = get_arg(args, "--gwps"))
out_dir   <- get_arg(args, "--out_dir", "genes")
max_beta  <- as.numeric(get_arg(args, "--max_beta", "-0.75"))
min_cells <- as.integer(get_arg(args, "--min_cells", "50"))
if (any(is.na(files)) || length(files) != 2) stop("need --essential and --gwps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ts <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

# ---- copied unchanged from inspre 1.0.1 ----
parse_hdf5_df <- function(hfile, entry = "obs") {
  colnames <- names(hfile[[entry]])
  cat_vals <- NULL
  if ("__categories" %in% colnames) {
    cat_names = names(hfile[[entry]][["__categories"]])
    get_cat <- function(cat_name) hfile[[entry]][["__categories"]][[cat_name]][]
    cat_vals = lapply(cat_names, get_cat)
    names(cat_vals) = cat_names
  }
  drop <- colnames %in% c("__categories")
  colnames <- colnames[!drop]
  get_col <- function(name) {
    if (name %in% names(cat_vals)) {
      indices = hfile[[entry]][[name]][]
      values = cat_vals[[name]][indices + 1]
    } else if ("categories" %in% names(hfile[[entry]][[name]])) {
      tryCatch({
        cat_vals[[name]] = hfile[[entry]][[name]][["categories"]][]
      }, warning = function(cond) print(name))
      indices = hfile[[entry]][[name]][["codes"]][]
      values = cat_vals[[name]][indices + 1]
    } else {
      values <- tryCatch(hfile[[entry]][[name]][],
                         error = function(cond) NA)
    }
    values
  }
  df <- as.data.frame(lapply(colnames, get_col))
  names(df) <- colnames
  df
}

calc_inst_effect_h5X <- function(inst_id, target, X, X_control, X_ids, X_vars) {
  inst_obs <- X_ids == inst_id
  n_target <- sum(inst_obs)
  n_control <- ncol(X_control)
  n_total <- n_target + n_control
  this_feature <- which(X_vars == target)
  X_exp_target <- X[this_feature, inst_obs]
  X_exp_control <- X_control[this_feature, ]
  Z <- c(rep(1, n_target), rep(0, n_control))
  inst_cor <- cor(Z, c(X_exp_target, X_exp_control))
  inst_beta <- cov(Z, c(X_exp_target, X_exp_control)) / var(Z)
  cor_se <- sqrt((1 - inst_cor^2) / (n_total - 2))
  as.data.frame(list(inst_id = inst_id, target = target, inst_beta = inst_beta,
                     inst_cor = inst_cor, cor_se = cor_se, n = n_total))
}
# ---- end of copied code ----

ntc <- "non-targeting"

select_screen <- function(screen, path) {
  ts(screen, ": opening ", path)
  h <- H5File$new(path, "r")
  on.exit(h$close_all())
  obs <- parse_hdf5_df(h, "obs")
  var <- parse_hdf5_df(h, "var")
  X <- h[["X"]]

  genes_guides <- filter(obs, gene != ntc) %>% distinct(gene_id, gene_transcript) %>%
    filter(gene_id %in% var$gene_id)
  cells_ntc <- obs$gene == ntc
  n_ntc <- sum(cells_ntc)
  ts(screen, ": ", nrow(obs), " cells, ", nrow(var), " genes, ", n_ntc,
     " non-targeting cells, ", nrow(genes_guides), " guides; loading controls")
  X_ntc <- X[, cells_ntc]
  rownames(X_ntc) <- var$gene_id

  ts(screen, ": guide effects")
  effects <- vector("list", nrow(genes_guides))
  for (k in seq_len(nrow(genes_guides))) {
    effects[[k]] <- calc_inst_effect_h5X(genes_guides$gene_transcript[k], genes_guides$gene_id[k],
                                         X, X_ntc, obs$gene_transcript, var$gene_id)
    if (k %% 500 == 0) ts(screen, ": ", k, "/", nrow(genes_guides), " guides")
  }
  guide_effects <- bind_rows(effects) %>%
    mutate(Z = inst_cor / cor_se, p = pt(Z, df = n - 2), p_adj = p.adjust(p, method = "fdr")) %>%
    arrange(target, p_adj) %>% filter(!duplicated(target)) %>%
    mutate(n_target_cells = n - n_ntc, screen = screen)
  write.csv(guide_effects, file.path(out_dir, sprintf("guide_effects_%s.csv", screen)), row.names = FALSE)

  kept <- filter(guide_effects, inst_beta < max_beta, n > n_ntc + min_cells)
  write.csv(kept, file.path(out_dir, sprintf("kept_%s.csv", screen)), row.names = FALSE)
  ts(screen, ": ", nrow(guide_effects), " targets, ", nrow(kept),
     " pass inst_beta < ", max_beta, " and > ", min_cells, " cells")
  list(kept = kept, var_order = var$gene_id)
}

res <- lapply(names(files), function(s) select_screen(s, files[[s]]))
names(res) <- names(files)

common <- intersect(res$essential$kept$target, res$gwps$kept$target)
common <- res$essential$var_order[res$essential$var_order %in% common]   # essential screen's gene order
write.csv(data.frame(common_genes = common), file.path(out_dir, "common_genes.csv"), row.names = FALSE)
ts("essential kept ", nrow(res$essential$kept), ", gwps kept ", nrow(res$gwps$kept),
   ", in both: ", length(common), " (the paper reports 521)")
