#!/usr/bin/env Rscript

# Checks that reading from .bin chunks gives exactly the same results as the original
# .Rdata path (get_info_from_reciprocal_plot_backup), and reports timings.
# Run from the repo root after tools/convert_tem_chunks.R:  Rscript tests/test_reciprocal_bin.R

suppressPackageStartupMessages({library(dplyr); library(Matrix)})
source("utils.R")

protein_index <- readRDS("data/tem_chunks/protein_index.rds")
PPI_sparse <- readRDS("data/PPI_sparse_sym.rds")
Biogrid_binary <- readRDS("data/Biogrid_binary.rds")
z_cache <- new.env()
rz_cache <- new.env()

n_proteins <- nrow(protein_index)
chunks <- as.integer(sub("\\.bin$", "", list.files("data/tem_chunks/z", pattern = "\\.bin$")))
chunks <- intersect(chunks, as.integer(sub("\\.bin$", "", list.files("data/tem_chunks/rz", pattern = "\\.bin$"))))
if (length(chunks) == 0) stop("No .bin chunks found. Run tools/convert_tem_chunks.R first.")

# 1. Every row/column of a few chunks matches the .Rdata chunk exactly
for (chunk_num in head(sort(chunks), 3)) {
  z <- readRDS(sprintf("data/tem_chunks/z/%d.Rdata", chunk_num))
  rz <- readRDS(sprintf("data/tem_chunks/rz/%d.Rdata", chunk_num))
  for (i in seq_len(nrow(z))) {
    stopifnot(identical(read_tem_block(sprintf("data/tem_chunks/z/%d.bin", chunk_num), i, n_proteins), as.numeric(z[i, ])))
    stopifnot(identical(read_tem_block(sprintf("data/tem_chunks/rz/%d.bin", chunk_num), i, n_proteins), as.numeric(rz[, i])))
  }
  cat(sprintf("chunk %d: all %d rows/columns identical\n", chunk_num, nrow(z)))
}

# 2. Full analysis output is identical to the original implementation
set.seed(1)
idx <- unlist(lapply(sample(chunks, min(5, length(chunks))), function(k) (k - 1) * 100 + sample(1:100, 2)))
idx <- idx[idx <= n_proteins]
t_new <- t_old <- 0
for (prot in protein_index$protein[idx]) {
  rm(list = ls(z_cache), envir = z_cache); rm(list = ls(rz_cache), envir = rz_cache)  # cold cache
  t_new <- t_new + system.time(new <- get_info_from_reciprocal_plot(prot))[["elapsed"]]
  t_old <- t_old + system.time(old <- get_info_from_reciprocal_plot_backup(prot))[["elapsed"]]
  stopifnot(identical(new, old))
}
cat(sprintf("%d proteins: output identical to original\n", length(idx)))
cat(sprintf("mean time per protein (cold): original %.3f s, .bin %.3f s\n", t_old / length(idx), t_new / length(idx)))
