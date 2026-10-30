#!/usr/bin/env Rscript

########################################################################################
###################### Convert xz-compressed tem chunks to uncompressed binary (.bin)
########################################################################################
# Each data/tem_chunks/{z,rz}/N.Rdata becomes N.bin next to it. The app reads a single
# protein's row/column from the .bin files with seek(), instead of decompressing the
# whole chunk. Chunks without a .bin file keep working through the .Rdata fallback.
#
# Layout: little-endian float64, one protein per contiguous block of n_proteins values.
#   z/N.bin  : rows of the z chunk     (block i = z_chunk[i, ])
#   rz/N.bin : columns of the rz chunk (block i = rz_chunk[, i])
#
# Run from the repo root:  Rscript tools/convert_tem_chunks.R

protein_index <- readRDS("data/tem_chunks/protein_index.rds")
n_proteins <- nrow(protein_index)

write_bin <- function(values, file) {
  con <- file(file, "wb")
  on.exit(close(con))
  writeBin(values, con, size = 8, endian = "little")
}

for (rdata in list.files("data/tem_chunks/z", pattern = "\\.Rdata$", full.names = TRUE)) {
  chunk_num <- sub("\\.Rdata$", "", basename(rdata))
  z_bin <- sprintf("data/tem_chunks/z/%s.bin", chunk_num)
  rz_rdata <- sprintf("data/tem_chunks/rz/%s.Rdata", chunk_num)
  rz_bin <- sprintf("data/tem_chunks/rz/%s.bin", chunk_num)

  if (!file.exists(rz_rdata)) {
    cat(sprintf("Skipping chunk %s: missing %s\n", chunk_num, rz_rdata))
    next
  }
  if (file.exists(z_bin) && file.exists(rz_bin)) next

  cat(sprintf("Converting chunk %s... ", chunk_num))
  z <- readRDS(rdata)
  rz <- readRDS(rz_rdata)
  stopifnot(ncol(z) == n_proteins, nrow(rz) == n_proteins)

  # R is column-major: t(z) puts each row of z in a contiguous block; rz columns already are
  write_bin(as.vector(t(z)), paste0(z_bin, ".tmp"))
  write_bin(as.vector(rz), paste0(rz_bin, ".tmp"))
  file.rename(paste0(z_bin, ".tmp"), z_bin)
  file.rename(paste0(rz_bin, ".tmp"), rz_bin)
  cat("done\n")
}
