#!/usr/bin/env Rscript

#**********************
# Command line to call the script in linux consol
#**********************

"Associate reads number to each corresponding region of peaks (number of reads from .bam files)

Usage:
  peaks_union_featureCounts.R [options] <rds_input> <bam1_input> <bam2_input>
  peaks_union_featureCounts.R -h | --help

Options:
  --output_rds <file>  Output rds file
  --output_txt <file>  Output matrix count txt file
  --output_csv <file> Output csv file
  -h, --help           Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)


#**********************
# Libraries and files loading
#**********************

suppressPackageStartupMessages({
  suppressWarnings(library(Rsubread))
  suppressWarnings(library(dplyr))
})
  

df = data.frame(readRDS(arguments$rds_input)) %>% 
      tibble::rownames_to_column(var = "number") %>%
      dplyr::mutate(GeneID = paste0("peak_", number), .keep = "unused") %>%
      dplyr::rename(Chr = seqnames, Start = start, End = end, Strand = strand) %>%
      dplyr::select(GeneID, Chr, Start, End, Strand) 

bam_files = c(arguments$bam1_input, arguments$bam2_input)
  
#**********************
# Script to create readcount_matrix
#**********************

# Count the number of reads, for each sample, present in the region resulting from the union
readCount <- featureCounts(files = bam_files,
                           annot.ext = df,  
                           isPairedEnd = TRUE,
                           nthreads = 1,
                           countChimericFragments = FALSE,
                           countMultiMappingReads = TRUE)
  

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
  count_df = tibble::rownames_to_column(data.frame(readCount$counts), "peakID")
  write.table(count_df, file = arguments$output_csv, sep = ";", row.names = FALSE) 
}



