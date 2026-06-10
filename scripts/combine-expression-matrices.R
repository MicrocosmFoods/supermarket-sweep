# Combine per-run expression matrices into a single matrix
#
# Each gzipped TSV in raw_data/rnaseq/ shares the same gene annotation
# columns (gene_id, gene_name, gene_biotype) and contributes its own set
# of per-sample columns (*_cpm, *_count). This script merges them on the
# annotation columns into one wide, gzipped matrix written to results/rnaseq/.

suppressPackageStartupMessages(library(data.table))

raw_dir     <- "raw_data/rnaseq"
results_dir <- "results/rnaseq"
key_cols    <- c("gene_id", "gene_name", "gene_biotype")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Runs to exclude from the combined matrix.
exclude <- c("44R9PF")

files <- list.files(raw_dir, pattern = "-expression-matrix\\.tsv\\.gz$", full.names = TRUE)
files <- files[!grepl(paste0("(", paste(exclude, collapse = "|"), ")-"), basename(files))]
stopifnot(length(files) > 0)
message("Combining ", length(files), " expression matrices:")
message(paste(" -", basename(files), collapse = "\n"))

tables <- lapply(files, function(f) fread(f, sep = "\t", header = TRUE))

# Merge all tables on the shared annotation columns.
combined <- Reduce(function(x, y) merge(x, y, by = key_cols, all = TRUE), tables)

out_file <- file.path(results_dir, paste0(Sys.Date(), "-combined-expression-matrix.tsv.gz"))
fwrite(combined, out_file, sep = "\t")

message("Wrote ", nrow(combined), " genes x ", ncol(combined), " columns to ", out_file)
