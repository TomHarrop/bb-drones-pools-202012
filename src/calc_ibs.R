#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(SNPRelate)
library(data.table)
library(ggplot2)

vcf_file <- snakemake@input[["vcf"]]
gds_file <- tempfile(fileext = ".gds")

# method = "biallelic.only": to exact bi-allelic and polymorhpic SNP data
# (excluding monomorphic variants); method = "copy.num.of.ref": to extract and
# store dosage (0, 1, 2) of the reference allele for all variant sites,
# including bi-allelic SNPs, multi-allelic SNPs, indels and structural variants.
snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only" )
gds_data <- snpgdsOpen(gds_file)

# SNPs on autosome with qual > 30
quals <- read.gdsn(index.gdsn(gds_data, "snp.annot/qual"))
autosomes <- startsWith(read.gdsn(index.gdsn(gds_data, "snp.chromosome")),
                        "NC_")

keep_snps <- which(quals > 30 & autosomes)

# run identity by state
ibs_res <- snpgdsIBS(gds_data,
                     snp.id = keep_snps,
                     autosome.only = FALSE)

# convert to df and write
ibs_df <- data.frame(ibs_res$ibs, row.names = ibs_res$sample.id)
colnames(ibs_df) <- ibs_res$sample.id

fwrite(ibs_df, 
       snakemake@output[["ibs"]],
       row.names = TRUE)

sessionInfo()
