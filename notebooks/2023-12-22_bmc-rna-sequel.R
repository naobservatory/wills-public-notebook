setwd("~/Dropbox/dev/work/playground/trial-tasks/")
library(tidyverse)
library(cowplot)
library(patchwork)
library(fastqcr)

# Compute approximate number of sequence bases from sequence length distribution table
compute_n_bases <- function(sld){
  min_lengths <- sub("-.*", "", sld$Length) %>% as.numeric
  max_lengths <- sub(".*-", "", sld$Length) %>% as.numeric
  mean_lengths <- (min_lengths + max_lengths)/2
  sld2 <- mutate(sld, Length_mean = mean_lengths)
  n_bases <- sum(sld2$Count * sld2$Length_mean)
  return(n_bases)
}

summarize_adapters <- function(adapter_tab){
  adapter_tab %>% select(-Position) %>% rowSums %>% max
}

# FASTQC processing functions
extract_qc_data_base <- function(qc){
  filename <- filter(qc$basic_statistics, Measure == "Filename")$Value
  n_reads <- filter(qc$basic_statistics, Measure == "Total Sequences")$Value %>% as.integer
  n_bases <- compute_n_bases(qc$sequence_length_distribution)
  pc_dup <- 100-qc$total_deduplicated_percentage
  per_base_quality <- filter(qc$summary, module == "Per base sequence quality")$status
  per_sequence_quality <- filter(qc$summary, module == "Per sequence quality scores")$status
  per_base_composition <- filter(qc$summary, module == "Per base sequence content")$status
  per_sequence_composition <- filter(qc$summary, module == "Per sequence GC content")$status
  per_base_n_content <- filter(qc$summary, module == "Per base N content")$status
  adapter_content <- filter(qc$summary, module == "Adapter Content")$status
  pc_adapters <- summarize_adapters(qc$adapter_content)
  tab_out <- tibble(n_reads = n_reads,
                    n_bases = n_bases,
                    pc_dup = pc_dup,
                    pc_adapters = pc_adapters,
                    per_base_quality = per_base_quality,
                    per_sequence_quality = per_sequence_quality,
                    per_base_composition = per_base_composition,
                    per_sequence_composition = per_sequence_composition,
                    per_base_n_content = per_base_n_content,
                    adapter_content = adapter_content,
                    filename = filename)
  return(tab_out)
}
extract_qc_data_paired <- function(qc, id_pattern, pair_pattern = ".*_(\\d)[_.].*"){
  base <- extract_qc_data_base(qc)
  filename <- base$filename[1]
  prefix <- sub(id_pattern, "", filename)
  read_pair <- sub(pair_pattern, "\\1", filename)
  tab_out <- base %>% mutate(sequencing_id = prefix, read_pair = read_pair) %>%
    select(sequencing_id, read_pair, everything())
  return(tab_out)
}

# Get fastqc filenames
data_dirs <- list(raw = "qc/raw", preproc = "qc/preproc", dedup = "qc/dedup", ribo = "qc/ribo")
id_patterns <- list(raw = "_masked_\\d.fastq.gz",
                   preproc = "_masked_fastp_\\d.fastq.gz",
                   dedup = "_masked_clumpify_\\d.fastq.gz",
                   ribo = "_masked_bbduk_\\d.fastq.gz")
qc_archives <- lapply(data_dirs, function(x) list.files(x, pattern = "*.zip", full.names = TRUE))
qc_data <- suppressMessages(lapply(qc_archives, function(x) lapply(x, function(y) qc_read(y, verbose = FALSE))))
qc_tab <- lapply(names(qc_data), function(x) 
  lapply(qc_data[[x]], function(y) extract_qc_data_paired(y, id_patterns[[x]])) %>% bind_rows %>% mutate(stage = x)) %>%
  bind_rows