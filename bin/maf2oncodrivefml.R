#!/usr/bin/env Rscript
library(dplyr)

# Convert a MAF to oncodrivefml input format

args <- commandArgs(trailingOnly = TRUE)

maf_file <- args[1]
output_file <- args[2]
exclude_file <- args[3]

if (is.na(maf_file) | is.na(output_file)) {
    stop("Usage: maf2oncodrivefml.R maf_file output_file")
} else if (! file.exists(maf_file)) {
    stop(paste("File", maf_file, "not found."))
}

maf <- read.table(maf_file, sep = "\t", header = T, quote = "", comment.char = "")
exclude_list <- if (!is.na(exclude_file)) read.table(exclude_file, header = F) else NULL

oncodrive_vars <- maf %>%
    select(Chromosome, POS_VCF, TUMOUR, REF_VEP, ALT_VEP) %>%
    rename("CHROMOSOME" = "Chromosome", "POSITION" = "POS_VCF", "SAMPLE" = "TUMOUR", "REF" = "REF_VEP", "ALT" = "ALT_VEP")

if (!is.null(exclude_list)) {
    oncodrive_vars <- oncodrive_vars %>%
    filter(!(SAMPLE %in% exclude_list$V1))
}

write.table(oncodrive_vars, file = output_file, row.names = F, sep = "\t", quote = F)
