
#########################
### Command line to call the script in linux consol
#########################

"Create plot for gene coverage

Usage:
  ATAC_gene_dynamics.R [options] <gr_dir> <readcount_dir>
  ATAC_gene_dynamics.R -h | --help

Options:
  -h, --help              Show this screen

" -> doc

library(docopt)
# arguments <- docopt(doc)

arguments <- docopt(doc, args=c(
  "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges",
  "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/macs2_output"
))



#########
### Initialization : libraries loading and function definitions
#########

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(GenomicRanges))
  suppressWarnings(library(stringr))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(ggthemes))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(gridExtra))
  suppressWarnings(library(tidyverse))
  suppressWarnings(library(rlang))  # pour utiliser case_when(!!!parse_exprs)
})

# # peaks_list_name = list.files(arguments$gr_dir, pattern = "_threshold_10_ann.gr.rds", full.names = TRUE)
# # # readcount_list_name = list.files(arguments$readcount_dir, pattern = "readcount", full.names = TRUE)

# 
# ########## utiliser la fonction prom_fun !!!! important !!!

readcount_name ="/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/macs2_output/Xvivo_00h_D1-D2-D3.readcount.csv"
peaks_name = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/Xvivo_00h_D1-D2-D3_threshold_10_ann.gr.rds"

peaks_reads_fun = function(peaks_name, readcount_name) {
  
  # Chargement des deux fichiers d'entrée, le fichier peaks a déjà subi un threshold
  readcount = read.table(readcount_name, sep = ";", header = TRUE)
  peaks = as.data.frame(readRDS(peaks_name)) %>%
    tibble::rownames_to_column(var = "sample_peak") %>%
    tidyr::separate(col = sample_peak, into = c("condition", "time", "donor", "peakID", "number"), sep = "_") %>%
    tidyr::unite(col = peakID, peakID, number, sep = "_", remove = TRUE) 
  
  # Fusionne les deux informations dans un Grange (liste de peaks et readcount)
  peaks = makeGRangesFromDataFrame(df = left_join(peaks, readcount, by = c("peakID", "condition", "time", "donor")) , 
                                   keep.extra.columns = TRUE, 
                                   ignore.strand = FALSE, 
                                   seqnames.field = "seqnames",
                                   start.field = "start", 
                                   end.field = "end", 
                                   strand.field = "strand")
  return(peaks)
}

# Sont donnés en entrée de code que les fichiers pour faire une cinétique
peaks_list_name = list("/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/MP_03h_D1-D2-D3_threshold_10_ann.gr.rds",
                       "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/MP_12h_D1-D2-D3_threshold_10_ann.gr.rds",
                       "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/MP_24h_D1-D2-D3-D5-D6_threshold_10_ann.gr.rds",
                       "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/Xvivo_00h_D1-D2-D3_threshold_10_ann.gr.rds")

readcount_list_name = list("/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/macs2_output/MP_03h_D1-D2-D3.readcount.csv",
                           "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/macs2_output/MP_12h_D1-D2-D3.readcount.csv",
                           "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/macs2_output/MP_24h_D1-D2-D3-D5-D6.readcount.csv",
                           "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/macs2_output/Xvivo_00h_D1-D2-D3.readcount.csv")


# ATAC - Création des granges annotés avec le nombre de reads par peaks (appel de fonction) #
#*******************************************************************************************#
peaks_reads_list = list()
for (i in 1:length(peaks_list_name)){
  name_sample = str_extract(peaks_list_name[[i]], "(?<=/genomic_ranges/).+(?=_threshold)")
  peaks_reads_list[[i]] = peaks_reads_fun(peaks_name = str_subset(peaks_list_name, name_sample),
                                         readcount_name = str_subset(readcount_list_name, name_sample))
  names(peaks_reads_list)[[i]] = name_sample
}

### Etude des peaks tombant dans l'espace intergénique ###
#********************************************************#

peaks_reads_intergenic = list()
for (i in 1:length(peaks_reads_list)) {
  peaks_reads_intergenic[[i]] = peaks_reads_list[[i]][which(peaks_reads_list[[i]]$Intergenic == TRUE)]
  names(peaks_reads_intergenic)[[i]] = names(peaks_reads_list)[[i]]
  # ajout de la colonne index 
  mcols(peaks_reads_intergenic[[i]])[ncol(mcols(peaks_reads_intergenic[[i]]))+1] = c(1:length(peaks_reads_intergenic[[i]]))
  colnames(mcols(peaks_reads_intergenic[[i]]))[ncol(mcols(peaks_reads_intergenic[[i]]))] = "index"
}

  # Union des timings : permet de lister toutes les régions où il y a au moins un pic itergénique à un moment donné
intergenic_gr_list_union = Reduce(GenomicRanges::union, peaks_reads_intergenic)
# Ajout de la colonne index
mcols(intergenic_gr_list_union) = c(1:length(intergenic_gr_list_union))
colnames(mcols(intergenic_gr_list_union)) = "index"

# Test de l'overlap des pics à chaque point de temps contre la liste générée avant.
# S'il y a une intersection des pics entre les différents points de temps, alors celle-ci a été fusionnée lors de l'union
# Il suffit alors de voir quels sont les pics qui ont participé à cette union, grâce à l'overlap des pics à chaque point de temps avec celle-ci

overlap_sum_up = tibble(union_ID = integer(),
                        union_ranges = character())

for (i in 1:length(peaks_reads_intergenic)){
  
  overlap = findOverlaps(intergenic_gr_list_union, peaks_reads_intergenic[[i]])
  
  link_df = tibble(union_ID = queryHits(overlap), new_col = subjectHits(overlap))
  colnames(link_df)[ncol(link_df)] = paste0("peakIndex_", names(peaks_reads_intergenic)[[i]])
  
  name_sample = names(peaks_reads_intergenic)[[i]]
  gr_df = data.frame(index = as.integer(peaks_reads_intergenic[[i]]$index),
                     peakID = peaks_reads_intergenic[[i]]$peakID,
                     nbreads = peaks_reads_intergenic[[i]]$nbreads)
  colnames(gr_df) = c(paste0("peakIndex_", name_sample),
                      paste0("peakID_", name_sample),
                      paste0("nbreads_", name_sample))

  union_df = data.frame(union_ID = as.integer(intergenic_gr_list_union$index),
                        start = as.data.frame(ranges(intergenic_gr_list_union))$start,
                        end = as.data.frame(ranges(intergenic_gr_list_union))$end) %>%
    tidyr::unite(col = "union_ranges", start, end, sep = "-", remove = TRUE)
    
  tab_recap = left_join(union_df, link_df, by = "union_ID")
  tab_recap = left_join(tab_recap, gr_df)

  overlap_sum_up = full_join(overlap_sum_up, tab_recap, by = c("union_ID", "union_ranges"))
}

overlap_sum_up[is.na(overlap_sum_up)] <- 0

# On remplace l'ID de chaque peak par sa valeur en nombre de reads 


## à mettre dans une fonction pour pouvoir faire la même chose pour chaque annotation ???





### Etude des peaks tombant dans les promoteurs FANTOM5 ###
#********************************************************#

peaks_reads_promoter = list()
for (i in 1:length(peaks_reads_list)) {
  peaks_reads_promoter[[i]] = peaks_reads_list[[i]][which(peaks_reads_list[[i]]$FANTOM5_promoter == TRUE)]
  names(peaks_reads_promoter)[[i]] = names(peaks_reads_list)[[i]]
}

# Union des timings : permet de lister toutes les régions où il y a au moins un pic itergénique à un moment donné
promoter_gr_list_union = Reduce(GenomicRanges::union, peaks_reads_promoter)
