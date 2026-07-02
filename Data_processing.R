
library(Matrix)
library(dplyr)


########################################################################################
###################### Data is already processed , so no need to run this code anymore
########################################################################################


########################################################################################
###################### Compressing size of row-z-normalized EffectiveConductance matrix file (tem.Rdata)
########################################################################################
{

# Load the matrix
tem <- readRDS("~/Downloads/FLOW/MASTERTRIANGLE/PPI/tem.Rdata")
cat(sprintf("Original dimensions: %d × %d\n", nrow(tem), ncol(tem)))
cat(sprintf("Original size (memory): %.2f GB\n", object.size(tem) / 1024^3))

# Get protein names
proteins <- rownames(tem)
n_proteins <- nrow(tem)

# Create directory structure
dir.create("data", showWarnings = FALSE)
dir.create("data/tem_chunks", showWarnings = FALSE)
dir.create("data/tem_chunks/z", showWarnings = FALSE)   # rows
dir.create("data/tem_chunks/rz", showWarnings = FALSE)  # columns

# Chunk parameters
chunk_size <- 100
n_chunks <- ceiling(n_proteins / chunk_size)

cat(sprintf("\nCreating %d chunks of ~%d proteins each...\n", n_chunks, chunk_size))

# Create chunks
for (chunk_num in 1:n_chunks) {
  
  # Calculate range
  start_idx <- (chunk_num - 1) * chunk_size + 1
  end_idx <- min(chunk_num * chunk_size, n_proteins)
  
  cat(sprintf("Processing chunk %02d: proteins %d-%d... ", 
              chunk_num, start_idx, end_idx))
  
  # Extract rows and columns for this chunk
  chunk_rows <- tem[start_idx:end_idx, , drop = FALSE]
  chunk_cols <- tem[, start_idx:end_idx, drop = FALSE]
  
  # Remove names to save space
  rownames(chunk_rows) <- NULL
  colnames(chunk_rows) <- NULL
  rownames(chunk_cols) <- NULL
  colnames(chunk_cols) <- NULL
  
  # Save rows (z)
  z_file <- file.path("data/tem_chunks/z", sprintf("%d.Rdata", chunk_num))
  saveRDS(chunk_rows, z_file, compress = "xz")
  
  # Save columns (rz)
  rz_file <- file.path("data/tem_chunks/rz", sprintf("%d.Rdata", chunk_num))
  saveRDS(chunk_cols, rz_file, compress = "xz")
  
  # Report sizes
  z_size_mb <- file.info(z_file)$size / 1024^2
  rz_size_mb <- file.info(rz_file)$size / 1024^2
  cat(sprintf("✓ z: %.1f MB, rz: %.1f MB\n", z_size_mb, rz_size_mb))
  
  # Clean up to save memory
  rm(chunk_rows, chunk_cols)
  gc()
}##Luka##

# Save protein index for quick lookup
protein_index <- data.frame(protein = proteins,index = 1:n_proteins,chunk = ceiling(1:n_proteins / chunk_size),stringsAsFactors = FALSE)

saveRDS(protein_index, "data/tem_chunks/protein_index.rds", compress = "xz")


}##Luka##


##################!!!!######################################################################
###################### Compressing size of row-z-normalized EffectiveConductance matrix file (tem.Rdata),with names
########################################################################################
{
  
  # Load the matrix
  tem <- readRDS("~/Downloads/FLOW/MASTERTRIANGLE/PPI/tem.Rdata")
  cat(sprintf("Original dimensions: %d × %d\n", nrow(tem), ncol(tem)))
  cat(sprintf("Original size (memory): %.2f GB\n", object.size(tem) / 1024^3))
  
  # Get protein names
  proteins <- rownames(tem)
  n_proteins <- nrow(tem)
  
  # Create directory structure
  dir.create("data", showWarnings = FALSE)
  dir.create("data/tem_chunks", showWarnings = FALSE)
  dir.create("data/tem_chunks/z", showWarnings = FALSE)   # rows
  dir.create("data/tem_chunks/rz", showWarnings = FALSE)  # columns
  
  # Chunk parameters
  chunk_size <- 100
  n_chunks <- ceiling(n_proteins / chunk_size)
  
  cat(sprintf("\nCreating %d chunks of ~%d proteins each...\n", n_chunks, chunk_size))
  
  # Create chunks
  for (chunk_num in 1:n_chunks) {
    
    # Calculate range
    start_idx <- (chunk_num - 1) * chunk_size + 1
    end_idx <- min(chunk_num * chunk_size, n_proteins)
    
    cat(sprintf("Processing chunk %02d: proteins %d-%d... ", 
                chunk_num, start_idx, end_idx))
    
    # Extract rows and columns for this chunk
    chunk_rows <- tem[start_idx:end_idx, , drop = FALSE]
    chunk_cols <- tem[, start_idx:end_idx, drop = FALSE]
    
    # Remove names to save space
    rownames(chunk_rows) <- NULL
    colnames(chunk_rows) <- NULL
    rownames(chunk_cols) <- NULL
    colnames(chunk_cols) <- NULL
    
    # Save rows (z)
    z_file <- file.path("data/tem_chunks/z", sprintf("%d.Rdata", chunk_num))
    saveRDS(chunk_rows, z_file, compress = "xz")
    
    # Save columns (rz)
    rz_file <- file.path("data/tem_chunks/rz", sprintf("%d.Rdata", chunk_num))
    saveRDS(chunk_cols, rz_file, compress = "xz")
    
    # Report sizes
    z_size_mb <- file.info(z_file)$size / 1024^2
    rz_size_mb <- file.info(rz_file)$size / 1024^2
    cat(sprintf("✓ z: %.1f MB, rz: %.1f MB\n", z_size_mb, rz_size_mb))
    
    # Clean up to save memory
    rm(chunk_rows, chunk_cols)
    gc()
  }##Luka##
  
  # Save protein index for quick lookup
  protein_index <- data.frame(protein = proteins,index = 1:n_proteins,chunk = ceiling(1:n_proteins / chunk_size),stringsAsFactors = FALSE)
  
  saveRDS(protein_index, "data/tem_chunks/protein_index.rds", compress = "xz")
  
  
}##Luka##



########################################################################################
###################### Compressing size of PPI score file (ECdz_min.Rdata)
########################################################################################
{
  thr <- 0.2
  PPI_matrix <- readRDS("~/Downloads/FLOW/MASTERTRIANGLE/PPI/ECdz_min.Rdata")
  lowest_score=min(PPI_matrix)
  
  PPI_matrix[PPI_matrix < thr] <- 0
  smat <- Matrix(PPI_matrix, sparse = TRUE)
  smat_sym <- forceSymmetric(smat, uplo = "U")
  saveRDS(smat_sym, "./data/PPI_sparse_sym.rds", compress = "xz")
  system("du -sh ./data/PPI_sparse_sym.rds")
  system("say done")
  
}##Luka##



########################################################################################
###################### Compressing size of CCI score file (CIM.Rdata)
########################################################################################
{
  thr <- 0
  CCI_matrix <- readRDS("~/Downloads/FLOW/MASTERTRIANGLE/PPI/CIM_CORUM+custommadepolymerase.Rdata")

  CCI_matrix[CCI_matrix < thr] <- 0

  smat <- Matrix(CCI_matrix, sparse = TRUE)
  smat_sym <- forceSymmetric(smat, uplo = "U")
  saveRDS(smat_sym, "./data/CCI_sparse_sym.rds", compress = "xz")
  system("du -sh ./data/CCI_sparse_sym.rds")
  system("say done")
  
}##Luka##

########################################################################################
###################### Compressing size of biogrid
########################################################################################
{
  thr <- 0
  matrix <- readRDS("~/Downloads/FLOW/MASTERTRIANGLE/PPI/interaction_matrix.Rdata")
  
  matrix <- matrix != 0
  matrix <- Matrix(matrix, sparse = TRUE)
  
  # Make symmetric by taking the union of upper and lower triangles
  # matrix <- matrix | t(matrix)
  matrix <- forceSymmetric(matrix, uplo = "U")
  
  saveRDS(matrix, "./data/Biogrid_binary.rds", compress = "xz")
}




