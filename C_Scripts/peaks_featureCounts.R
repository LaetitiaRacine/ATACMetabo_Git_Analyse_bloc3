#!/usr/bin/env Rscript

#**********************
# Command line to call the script in linux consol
#**********************

"Associate reads number to each corresponding region of peaks (number of reads from .bam files)

Usage:
  peaks_featureCounts.R [options] <csv_input> <bam_input>...
  peaks_featureCounts.R -h | --help

Options:
  --output_rds <file>  Output rds file
  --output_txt <file>  Output matrix count txt file
  --output_csv <file>  Output csv file
  -h, --help           Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)


#**********************
# Libraries loading
#**********************

suppressPackageStartupMessages({
  suppressWarnings(library(Rsubread))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(stringr))
})

df = read.table(arguments$csv_input, sep=";", header = TRUE, stringsAsFactors = FALSE)

if (str_detect(arguments$csv_input, pattern = "_ann")) {
  df = df[,c(1:4,6)] %>%
    dplyr::mutate(Name = str_extract(string = arguments$csv_input, pattern = "(?<=/genomic_ranges/).+(?=_ann)")) %>%
    dplyr::mutate(GeneID = paste(Name, X, sep="_peak_") ) %>%
    dplyr::rename(Chr = seqnames, Start = start, End = end, Strand = strand) %>%
    dplyr::select(-Name, - X) %>%
    dplyr::select(GeneID, Chr, Start, End, Strand) 
}

#**********************
# Script to create readcount_matrix
#**********************
  
readCount <- featureCounts(files = arguments$bam_input,
                           annot.ext = df,  # Besoin d'avoir les noms GeneID, Chr, Start, End et Strand dans le tableau pour fonctionner
                           isPairedEnd = TRUE,
                           nthreads = 1,
                           countChimericFragments = FALSE,
                           countMultiMappingReads = TRUE
)

if (!is.null(arguments$output_rds)) {
  print("Save rds...")
  saveRDS(readCount, file = arguments$output_rds)
}

if (!is.null(arguments$output_txt)) {
  print("Save txt...")
  matrix_count = readCount$counts
  write.table(matrix_count, file = arguments$output_txt, sep = "\t", quote = FALSE)
}

if (!is.null(arguments$output_csv)) {
  print("Save csv...")
  count_df = tibble::rownames_to_column(data.frame(readCount$counts), "name")
  count_df = tidyr::separate(count_df, col = name, into = c("condition","time","donor", "peak", "num"), sep="_", remove = TRUE) 
  count_df = tidyr::unite(count_df, col = peakID, peak, num, sep="_", remove = TRUE)
  if (!str_detect(arguments$csv_input, pattern = "_ann")) { colnames(count_df) = c("condition","time","donor","peakID","nbreads") }
  write.table(count_df, file = arguments$output_csv, sep = ";", row.names = FALSE) 
  }
