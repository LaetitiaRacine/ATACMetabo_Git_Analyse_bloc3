#!/usr/bin/env Rscript

#**********************
# Command line to call the script in linux consol
#**********************

"Create and manipulate GRanges files

Usage:
  GRanges.R [options] from_csv <csv_file>
  GRanges.R [options] union <gr_file> <gr_file>...
  GRanges.R [options] intersect <gr_file> <gr_file>...
  GRanges.R -h | --help

Options:
  -o, --output <gr_file>  Output file
  -h, --help              Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#**********************
# Libraries loading
#**********************

suppressPackageStartupMessages({
  suppressWarnings(library(GenomicRanges))
  suppressWarnings(library(dplyr))
})

#**********************
# One code per option
#**********************


if (arguments$from_csv) { 
  ### Transform df.csv file from broadPeak into GRange object
  df_csv = read.table(arguments$csv_file, sep=";", header = TRUE, stringsAsFactors = FALSE)
  df_csv$Chr = as.factor(df_csv$Chr)
  chromosomes_represented = levels(as.factor(df_csv$Chr))
  hg19_seqlengths = c("chr1"=249250621, "chr10"=135534747, "chr11"=135006516, "chr12"=133851895, "chr13"=115169878,
                      "chr14"=107349540, "chr15"=102531392, "chr16"=90354753, "chr17"=81195210, "chr18"=78077248,
                      "chr19"=59128983, "chr2"=243199373, "chr20"=63025520, "chr21"=48129895, "chr22"=51304566,
                      "chr3"=198022430, "chr4"=191154276, "chr5"=180915260, "chr6"=171115067, "chr7"=159138663,
                      "chr8"=146364022, "chr9"=141213431, "chrMT"=16571, "chrX"=155270560, "chrY"=59373566)
  hg19_seqlengths = hg19_seqlengths[names(hg19_seqlengths) %in% chromosomes_represented]
  gr = GRanges(seqnames = df_csv$Chr,
                ranges = IRanges(df_csv$Start, df_csv$End),
                strand = "*")
  names(gr) = df_csv$GeneID # add peak names to Grange (sample_peak_x)
  genome(gr) = "hg19" # add genome information to the Grange
  seqlengths(gr) = hg19_seqlengths # add chr length to each chr
  saveRDS(gr, file=arguments$output)
}


if (arguments$intersect) {
  # Determine common peaks in a list of Granges (intersection)
  gr_list = lapply(X = arguments$gr_file, FUN = readRDS)
  gr_intersection = Reduce(GenomicRanges::intersect, gr_list)
  saveRDS(gr_intersection, file=arguments$output)
}

if (arguments$union) {
  # Determine all peaks in a list of Granges (union)
  gr_list = lapply(arguments$gr_file, readRDS)
  gr_intersection = Reduce(GenomicRanges::union, gr_list)
  saveRDS(gr_intersection, file=arguments$output)
}
