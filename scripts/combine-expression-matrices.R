# Combine per-run expression matrices into a single matrix
#
# Each TSV in raw_data/ shares the same gene annotation columns
# (gene_id, gene_name, gene_biotype) and contributes its own set of
# per-sample columns (*_cpm, *_count). This script merges them on the
# annotation columns into one wide matrix written to results/.

suppressPackageStartupMessages(library(data.table))

raw_dir     <- "raw_data"
results_dir <- "results"
key_cols    <- c("gene_id", "gene_name", "gene_biotype")

files <- list.files(raw_dir, pattern = "-expression-matrix\\.tsv$", full.names = TRUE)
stopifnot(length(files) > 0)
message("Combining ", length(files), " expression matrices:")
message(paste(" -", basename(files), collapse = "\n"))

tables <- lapply(files, function(f) fread(f, sep = "\t", header = TRUE))

# Merge all tables on the shared annotation columns.
combined <- Reduce(function(x, y) merge(x, y, by = key_cols, all = TRUE), tables)

out_file <- file.path(results_dir, paste0(Sys.Date(), "-combined-expression-matrix.tsv"))
fwrite(combined, out_file, sep = "\t")

message("Wrote ", nrow(combined), " genes x ", ncol(combined), " columns to ", out_file)
