#**********************
# Command line to call the script in linux consol
#**********************

"Add annotations in GRanges files

Usage:
  annotate_grange.R [--output_csv <csv_file>] -a <annotations_file> <grange_input_file> <grange_annotated_file>
  annotate_grange.R -h | --help

Options:
  --output_csv <csv_file>  Output annotated file
  -h, --help               Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#**********************
# Libraries loading and function definition
#**********************

suppressPackageStartupMessages({
  suppressWarnings(library(GenomicRanges))
  suppressWarnings(library(dplyr))
})

loadRData<-function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls()!="fileName"])
}

#**********************
# Files loading
#**********************

# the "big" file containing "all" the annotations
all_annotations = loadRData(arguments$annotations_file)

# the grange to annotate
gr = readRDS(arguments$grange_input_file)

#**********************
# Annotation of the grange
#**********************

# First a matrix is created filled with FALSE and added to the Grange
annotations_types = levels(factor(all_annotations$annotation))
metadata = matrix(FALSE, ncol = length(annotations_types), nrow = length(gr))
colnames(metadata) = annotations_types
mcols(gr) = metadata # ajout de metadata : une colonne par annotation toutes à false au démarrage

# for each of the annotations types an overlap is calculated and used to assigned the peak as TRUE when overlapping with the annotation
for (i in 1:ncol(metadata)){
  sub_annot = all_annotations[all_annotations$annotation == annotations_types[i]]
  overlaps = findOverlaps(gr, sub_annot)
  mcols(gr)[queryHits(overlaps),i] = TRUE
}

colnames(mcols(gr)) = case_when(colnames(mcols(gr)) == "3' UTR" ~ "UTR3P",
                                colnames(mcols(gr)) == "5' UTR" ~ "UTR5P",
                                colnames(mcols(gr)) == "CpG Island" ~ "CpG",
                                colnames(mcols(gr)) == "CTCF" ~ "CTCF",
                                colnames(mcols(gr)) == "EXON" ~ "Exons",
                                colnames(mcols(gr)) == "FANTOM5_promoter" ~ "FANTOM5_promoter",
                                colnames(mcols(gr)) == "INTRON" ~ "Introns",
                                colnames(mcols(gr)) == "Promoter_+-1000"  ~ "TSS_mp1kb")

mcols(gr) = as_tibble(mcols(gr)) %>%
  dplyr::mutate(Intergenic = ifelse(UTR3P == FALSE & UTR5P == FALSE & Exons == FALSE & Introns == FALSE & FANTOM5_promoter == FALSE & TSS_mp1kb == FALSE, TRUE, FALSE)) %>%
  dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_Intergenic = ifelse(Intergenic == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_in_intron = ifelse(Introns == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_in_exon = ifelse(Exons == TRUE & CTCF == TRUE, TRUE, FALSE)) 

saveRDS(gr, arguments$grange_annotated_file)
if (!is.null(arguments$output_csv)) {
  write.csv2(as.data.frame(gr), arguments$output_csv)
}
