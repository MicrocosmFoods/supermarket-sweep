library(tidyverse)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(BiocParallel)

###############################
# RNA-seq analysis of full supermarket sweep data
###############################

# sample map
sample_map <- read_tsv("metadata/rnaseq/2026-09-25-rnaseq-sample-map.tsv", col_names = TRUE) %>% 
  set_names(c("sample_code", "batch", "sample_id", "lps_addition", "pre_treatment", "sample_type", "main_ingredient", "food_type", "replicate")) %>% 
  mutate(sample = paste0(sample_id, "_count")) %>% 
  select(sample, batch, sample_type, replicate)

# counts
sample_counts <- read_tsv("results/rnaseq/2026-09-25-combined-expression-matrix.tsv") %>% 
  filter(gene_biotype == "protein_coding") %>% 
  select(gene_id, gene_name, ends_with("_count"))

# check sample map and sample counts that they match
map_samples <- sample_map$sample
count_samples <- setdiff(colnames(sample_counts), c("gene_id", "gene_name"))
# in sample map but not the count matrix
setdiff(map_samples, count_samples) %>% sort()
# in count matrix but not sample map
setdiff(count_samples, map_samples) %>% sort()

# remove singleton samples that have no additional replicates
singleton_samples <- sample_map %>% 
  distinct() %>% 
  group_by(sample_type) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  pull(sample)

# filtered sample counts and turn into a matrix
filtered_sample_counts_matrix <- sample_counts %>% 
  select(-any_of(singleton_samples)) %>% 
  pivot_longer(
    !c(gene_id, gene_name),
       names_to = "sample",
       values_to = "count") %>% 
  mutate(count = round(count)) %>% 
  select(gene_id, sample, count) %>% 
  pivot_wider(names_from = sample, values_from = count) %>% 
  column_to_rownames("gene_id") %>% 
  as.matrix()

sample_metadata <- sample_map %>% 
  distinct(sample, sample_type, batch, replicate) %>% 
  dplyr::slice(match(colnames(filtered_sample_counts_matrix), sample)) %>% 
  column_to_rownames("sample") %>% 
  mutate(sample_type = factor(make.names(sample_type)),
         replicate = factor(replicate),
         batch = factor(batch))

# check metadata and counts matrix match
all(rownames(sample_metadata) == colnames(filtered_sample_counts_matrix))
levels(sample_metadata$sample_type)

# positive control LPS+ as reference
sample_metadata$sample_type <- relevel(sample_metadata$sample_type,
                                       ref="positive.control")

# DESeq object
dds <- DESeqDataSetFromMatrix(filtered_sample_counts_matrix,
                              colData = sample_metadata,
                              design = ~ sample_type)

# remove low counts
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]

# use BiocParallel for constructing dds object
register(MulticoreParam(workers = parallel::detectCores() - 2))
dds <- DESeq(dds, parallel = TRUE)

# QC checks
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("batch"), ntop=500)
plotPCA(vsd, intgroup = "sample_type")



pd <- plotPCA(vsd, intgroup = c("sample_type", "batch"), returnData = TRUE)
ggplot(pd, aes(PC1, PC2, color = sample_type == "negative.control", shape = batch)) + geom_point(size = 2) + scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 4))

ggplot(pd, aes(PC1, PC2)) +
  geom_point(color = "grey80") +
  geom_point(data = filter(pd, sample_type %in% c("positive.control", "negative.control")),
             aes(color = batch, shape = sample_type), size = 3)

