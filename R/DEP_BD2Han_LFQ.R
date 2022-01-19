#! R

##script to run on LFQ data provided herein
##NB that data was parsed first using parse_BD2Han_LF!_141220.sh

##based on DEP: https://www.bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html

library("DEP")
library("tidyverse")
library("limma")
library("SummarizedExperiment")

source("R/DEP_BD2Han_LFQ.func.R")
data_1 <- as.data.frame(readr::read_tsv("data/LFQ/BD2Han_LFQ_141220.parsed.tsv"))
data_2 <- as.data.frame(readr::read_tsv("data/LFQ/Endometrial_Proteins_Bruce_220121.parsed.tsv"))

sample_map <- readr::read_tsv("data/metadata/sample_map.tsv") %>%
              dplyr::mutate(group_base = unlist(lapply(group_ID, function(f){
                paste(strsplit(unlist(f), "_")[[1]][1:2], collapse="_")
                })))

##remove NA in UniprotID
data_clean <- dplyr::filter(.data = data_2, !is.na(UniprotID)) %>%
              dplyr::filter(!is.na(LFQ_intensity_R9_2))

##those cols with uniprot but not Gene_names, use UniprotID in Gene_names
##these aren't found in Biomart either!
data_clean_ng <- dplyr::mutate(.data = data_clean, Gene_names = if_else(is.na(Gene_names), UniprotID, Gene_names))

multis <- test_multi(data_clean_ng, "Gene_names") %>%
          dplyr::select(Gene_names) %>% unlist()

##remove those multi Gene_names, NA for Gene_names
data_munique <- dplyr::filter(.data = data_clean_ng, Gene_names %in% multis) %>%
                dplyr::mutate(Gene_names = paste0(Gene_names, "_", UniprotID))
data_unique <- dplyr::filter(.data = data_clean_ng, !Gene_names %in% multis)
data_uniq <- rbind(data_munique, data_unique)
tm <- test_multi(data_uniq, "Gene_names")
dim(tm)

##if multis remain...
if(dim(tm) > 0){
  mm <- max_multi(data_uniq, "Gene_names")
  data_uniq_1 <- data_uniq %>% dplyr::filter(! Gene_names %in% unique(tm$Gene_names)) %>% as.data.frame()
  data_uniq <- rbind(data_uniq_1, mm)
}
test_multi(data_uniq, "Gene_names")

##make unique
data_uni <- make_unique(data_uniq, "Gene_names", "UniprotID", delim = ";")

##make_unique
colnames(data_uni) <- gsub("LFQ_intensity_", "", colnames(data_uni))
colnames(data_uni) <- gsub("R", "R_", colnames(data_uni))
colnames(data_uni) <- gsub("__", "", colnames(data_uni))

##save that
data_uni_tb <- data_uni %>% dplyr::select(name, ID, dplyr::everything()) %>% as_tibble()
dir.create("output/LFQ", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(data_uni_tb, file = "output/LFQ/BD2Han.protein_LFQ_220121.parse_clean.csv")

##make SE
LFQ_columns <- grep("R_", colnames(data_uni))
raw <- data_uni[, LFQ_columns]
raw[raw == 0] <- NA
raw <- as.matrix(log2(raw))
rownames(raw) <- data_uni$name

proteins_unique <- data_uni
rownames(proteins_unique) <- proteins_unique$name

row_data <- proteins_unique[, -LFQ_columns]
row_data <- row_data[,-c(1,2,3)]
rownames(row_data) <- row_data$name

run_expt_group(sample_map, group = "group_ID", raw, row_data, outdir = "output/LFQ")
run_expt_group(sample_map, group = "group_base", raw, row_data, outdir = "output/LFQ")
run_expt_group(sample_map, group = "group_ID", raw, row_data, exclude_group = c("T_STD_HI", "T_LG"), outdir = "output/LFQ")
run_expt_group(sample_map, group = "group_ID", raw, row_data, convert_group = c("T_LG","T_STD_HI"), outdir = "output/LFQ")
