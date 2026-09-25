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

# Samples to drop (both their _cpm and _count columns are removed)
exclude_samples <- c("44R9PF_17", "44R9PF_18", "44R9PF_19", "44R9PF_20", "44R9PF_21")

files <- list.files(raw_dir, pattern = "-expression-matrix\\.tsv\\.gz$", full.names = TRUE)
stopifnot(length(files) > 0)
message("Combining ", length(files), " expression matrices:")
message(paste(" -", basename(files), collapse = "\n"))

tables <- lapply(files, function(f) fread(f, sep = "\t", header = TRUE))

# Merge all tables on the shared annotation columns.
combined <- Reduce(function(x, y) merge(x, y, by = key_cols, all = TRUE), tables)

# Drop excluded samples
drop_cols <- as.vector(outer(exclude_samples, c("_cpm", "_count"), paste0))
missing   <- setdiff(drop_cols, names(combined))
if (length(missing) > 0) warning("Columns not found, so not dropped: ", paste(missing, collapse = ", "))
drop_cols <- intersect(drop_cols, names(combined))
if (length(drop_cols) > 0) combined[, (drop_cols) := NULL]
message("Dropped ", length(drop_cols), " columns for ", length(exclude_samples), " excluded samples")

out_file <- file.path(results_dir, paste0(Sys.Date(), "-combined-expression-matrix.tsv.gz"))
fwrite(combined, out_file, sep = "\t")

message("Wrote ", nrow(combined), " genes x ", ncol(combined), " columns to ", out_file)
