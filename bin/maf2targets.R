#!/usr/bin/env Rscript

# Split a DERMATLAS keep MAF into per-sample region files for bcftools --targets-file.
# Writes sbs-dbs_targets_<sample>.tsv and id_targets_<sample>.tsv into the current
# directory and a samples.tsv manifest listing unique Tumor_Sample_Barcode values.

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
})

write_tsv <- function(df, sample, type) {
    write.table(df, paste0(type, "_", sample, ".tsv"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
maf_file <- args[1]

if (is.na(maf_file)) {
    stop("Usage: maf2targets.R <keep.maf>")
}

maf <- read.table(maf_file, header = TRUE, sep = "\t", check.names = FALSE)

chr_order <- str_sort(unique(maf$Chromosome), numeric = TRUE)
samples <- unique(maf$Tumor_Sample_Barcode)

writeLines(samples, "samples.tsv")

for (sample in samples) {
    maf_vars <- maf %>%
        filter(Tumor_Sample_Barcode == sample) %>%
        select(Chromosome, POS_VCF, End_Position, Variant_Type) %>%
        distinct() %>%
        arrange(factor(Chromosome, levels = chr_order))
    sbs_dbs <- maf_vars %>% filter(Variant_Type %in% c("SNP", "DNP"))
    indel <- maf_vars %>% filter(Variant_Type %in% c("DEL", "INS"))
    write_tsv(sbs_dbs, sample, "sbs-dbs_targets")
    write_tsv(indel, sample, "id_targets")
}
