#! R

##script to run on LFQ data provided herein
##NB that data was parsed first using parse_BD2Han_LF!_141220.sh

##based on DEP: https://www.bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html

library("DEP")
library("tidyverse")
library("limma")
library("SummarizedExperiment")

p1 <- readr::read_csv("../data/LFQ/Endometrial PhosphoProteins_Bruce_220121.csv")
sample_map <- readr::read_tsv("../data/metadata/sample_map.tsv") %>%
              dplyr::mutate(group_base = unlist(lapply(group_ID, function(f){
                paste(strsplit(unlist(f), "_")[[1]][1:2], collapse="_")
                })))

gns <- unlist(lapply(p1$`Fasta headers`, function(f){strsplit(strsplit(f,"GN=")[[1]][2], "[ ;]")[[1]][1]}))
trfhs <- unlist(lapply(p1$`Fasta headers`, function(f){strsplit(f,"\\|")[[1]][2]}))

##add new cols
p1$Gene_names <- paste0(gns, "_", trfhs)
p1$Uniprot_ID <- trfhs

##colnames horrorshow
outdir <- "output/LFQ"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
colnames(p1) <- gsub("LFQ intensity ", "R_", colnames(p1))
p1_tb <- p1 %>% dplyr::select(Gene_names, Uniprot_ID, dplyr::everything())
readr:::write_csv(p1_tb, file = paste0(outdir,"/BD2Han.phosphoprotein_220121.parse.csv"))
