
# DEseq : on étudie les peaks présents dans les deux conditions testés
# Les peaks sont-ils plus grands ou plus petits ? 

#**************************
# Libraries and directory
#**************************

dir = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/"
dir_output = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/D_Analysis/ATAC_differential_accessibility/"

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls()!="fileName"])
}

library(DESeq2)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(tibble)
library(tidyr)
library(ggplot2)


#####################################################################################
################# CODE MANUEL POUR UN SEUL ECHANTILLON ##############################  => à adapter pour snakefile
#####################################################################################


#**************************************
#* Chargement des fichiers nécessaires
#**************************************

# Pour faire l'analyse différentielle entre deux conditions, on a besoin : 
# - le Grange de l'union des conditions (donneurs fusionnés avec threshold pour chaque condition)
# gr_union = readRDS(paste0(dir_output, "2DG_03h_D1-D2-D3_vs_MP_03h_D1-D2-D3_threshold_10_ann.gr.rds"))
name_gr_union = paste0(dir_output, "2DG_24h_D1-D2-D3_vs_MP_24h_D1-D2-D3-D5-D6_threshold_10_ann.gr.rds")
gr_union = readRDS(name_gr_union)

cond1 = str_extract(name_gr_union, pattern = "(?<=_vs_)[:alnum:]{2,5}_[:digit:]{2}h(?=_D)")
cond2 = str_extract(name_gr_union, pattern = "(?<=differential_accessibility/)[:alnum:]{2,5}_[:digit:]{2}h(?=_D)")
name = paste0(cond2, "_vs_", cond1)

# - le bam de chaque condition pour tous les donneurs (avant passage au bloc2, un fichier par donneur)
name_bam_cond1 = str_subset(list.files(path = "A_Initial_data/merged_analysis/bam", pattern = cond1, full.names = TRUE), 
                            pattern = "_D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
name_bam_cond2 = str_subset(list.files(path = "A_Initial_data/merged_analysis/bam", pattern = cond2, full.names = TRUE), 
                            pattern = "_D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
bam_files = c(name_bam_cond1, name_bam_cond2)

# - le Grange de chaque condition en donneurs fusionnés avec threshold
gr_cond1 = readRDS(str_subset(list.files(path = "A_Initial_data/merged_analysis/genomic_ranges", pattern = cond1, full.names = TRUE), 
                                   pattern = "_threshold_[:digit:]{2,3}_ann.gr.rds"))
gr_cond2 = readRDS(str_subset(list.files(path = "A_Initial_data/merged_analysis/genomic_ranges", pattern = cond2, full.names = TRUE), 
                                   pattern = "_threshold_[:digit:]{2,3}_ann.gr.rds"))

#**********************
# Create readcount matrix
#**********************

# Mise en forme du Grange en data.frame pour le readcount
df = data.frame(gr_union) %>%
  tibble::rownames_to_column(var = "number") %>%
  dplyr::mutate(GeneID = paste0("peak_", number), .keep = "unused") %>%
  dplyr::rename(Chr = seqnames, Start = start, End = end, Strand = strand) %>%
  dplyr::select(GeneID, Chr, Start, End, Strand)

# Calcul du readcount
readCount <- featureCounts(files = bam_files,
                           annot.ext = df,  
                           isPairedEnd = TRUE,
                           nthreads = 1,
                           countChimericFragments = FALSE,
                           countMultiMappingReads = TRUE)

df_annot = readCount$annotation

# Enregistrement des readcount
print("Save txt...")
matrix_count = readCount$counts
write.table(matrix_count, 
            file = paste0(dir_output, "test_featurecounts_", name, ".txt"),
            sep = "\t", quote = FALSE)

print("Save csv...")
count_df = tibble::rownames_to_column(data.frame(readCount$counts), "GeneID")
write.table(count_df, 
            file = paste0(dir_output, "test_readcount_", name, ".csv"),
            sep = ";", row.names = FALSE)

df_readcount = left_join(df_annot, count_df, by = "GeneID") %>%
  dplyr::mutate(region = str_replace(GeneID, "peak", "region"), .before = Chr) %>%
  dplyr::select(-GeneID) 
  
#**********************
# Look on the readcount matrix
#**********************

# Ajout de l'information peak ou pas peak
# Lors du readcount, on perd l'info de la présence d'un peak détecté par peak calling ou non 
# On va faire des overlap pour voir si les régions extraites de l'union correspondent à un peak dans chaque condition
union_gr = makeGRangesFromDataFrame(df = df_readcount, 
                                   keep.extra.columns = TRUE,
                                   start.field = "Start",
                                   end.field = "End",
                                   seqnames.field = "Chr",
                                   ignore.strand = FALSE)

tab = data.frame(peak_cond1 = as.logical(countOverlaps(union_gr, gr_cond1)),
                 peak_cond2 = as.logical(countOverlaps(union_gr, gr_cond2))) %>%
  tibble::rownames_to_column(var = "region") %>%
  dplyr::mutate(region = paste0("region_", region))

df_readcount = left_join(df_readcount, tab, by = "region")


# Pour les volcano, on garde seulement les régions qui correspondent à un peak dans les deux conditions 

df_peaks_persistant = df_readcount %>% 
  dplyr::filter(peak_cond1 == TRUE & peak_cond2 == TRUE) %>%
  dplyr::select(-peak_cond1, -peak_cond2)

df_type_peak = data.frame(differential = name,
                          persistant_peaks = nrow(df_peaks_persistant),
                          close_open_peaks = nrow(df_readcount %>% dplyr::filter(peak_cond1 == FALSE & peak_cond2 == TRUE)),
                          open_close_peaks = nrow(df_readcount %>% dplyr::filter(peak_cond1 == TRUE & peak_cond2 == FALSE)))

matrix_peaks_volcano = df_peaks_persistant[,7:ncol(df_peaks_persistant)]

#**********************
# Deseq on readcount matrix
#**********************

# DEseq2 parameters (chaque sample a 3 donneurs donc 6 lignes au total)
coldata <- data.frame(sample = str_extract(bam_files, pattern = "(?<=bam/)[:alnum:]{2,5}_[:digit:]{2}h_D[:digit:]{1}(?=.bam)")) %>%
  dplyr::mutate(condition = ifelse(str_detect(sample, pattern = cond1), "before", "after"), 
                type = "paired-end")

# Run DEseq
# Dans coldata, une ligne correspond à une colonne de count_df
# Normalisation des données avec DEseq
dds <- DESeqDataSetFromMatrix(
  countData = matrix_peaks_volcano,
  colData = coldata,
  design = ~ condition)

## Setting the samples tagged "before" as reference
dds$condition <- relevel(dds$condition, ref = "before") 
dds <- DESeq(dds)

# Affichage de la matrice normalisée
cm = data.frame(counts(dds, normalized=TRUE))

res <- results(dds)
res = as_tibble(res)
res = res %>% 
  mutate(region = paste0(name,"_region_", 1:nrow(res))) %>%
  relocate(region, .before = baseMean)

#**********************
# Annotate DESeq2 results
#**********************

# 1 Annotate the region resulting from the union of peaks at two consecutive times ## Idem script bloc 2 : annotate_grange.R

all_annotations = loadRData(paste0(dir, "A_Initial_data/Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_woThisto_FANTOM5_prom_gr.rda"))
FANTOM_prom_gr = loadRData(paste0(dir, "A_Initial_data/prom_gene_fantom_gr.rdata"))
annotations_types = levels(factor(all_annotations$annotation))

# First a matrix is created filled with FALSE and added to the Grange
metadata = matrix(FALSE, ncol = length(annotations_types), nrow = length(gr_union))
colnames(metadata) = annotations_types
mcols(gr_union) = metadata

# for each of the annotations types an overlap is calculated and used to assigned the peak as TRUE when overlapping with the annotation
for (i in 1:ncol(metadata)){
  sub_annot = all_annotations[all_annotations$annotation == annotations_types[i]]
  overlaps = findOverlaps(gr_union, sub_annot)
  mcols(gr_union)[queryHits(overlaps),i] = TRUE
}

# Adding overlap with fantom5 promoter database
overlap_with_FANTOM = findOverlaps(gr_union, FANTOM_prom_gr)
mcols(gr_union)[,"FANTOM5_promoter"] = FALSE
mcols(gr_union)[queryHits(overlap_with_FANTOM), "FANTOM5_promoter"] = TRUE

colnames(mcols(gr_union)) = c("UTR3P","UTR5P","CpG", "CTCF","Exons","Introns","TSS_mp1kb","FANTOM5_promoter")

mcols(gr_union) = as_tibble(mcols(gr_union)) %>%
  dplyr::mutate(Intergenic = ifelse(UTR3P == FALSE & UTR5P == FALSE & Exons == FALSE & Introns == FALSE & FANTOM5_promoter == FALSE & TSS_mp1kb == FALSE, TRUE, FALSE)) %>%
  dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_Intergenic = ifelse(Intergenic == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_in_intron = ifelse(Introns == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
  dplyr::mutate(CTCF_in_exon = ifelse(Exons == TRUE & CTCF == TRUE, TRUE, FALSE))

df = as_tibble(gr_union) %>% 
  tibble::rownames_to_column(var = "number") %>%
  dplyr::mutate(region = paste0(name, "_region_", number)) %>%
  dplyr::select(-number) %>%
  relocate(region, .before = seqnames)

# 2 Adding the information concerning region resulting from the union to DEseq results

DEseq_results = res %>% 
  mutate(regulation = case_when(
      pvalue < 0.01 ~ "significative",
      pvalue > 0.01 ~ "Non-significative"))

# Adding start|end|seqnames informations to regions
DEseq_results_annotated = left_join(DEseq_results,  df, by = "region") 

## enlever les NA ! 

write.csv2(x = DEseq_results_annotated, file = paste0(dir_output, "/DEseq_results_", name, ".csv"))

df_recap_peak = table(DEseq_results_annotated$regulation)   
  

#**********************
# Volcano_plot
#**********************

DEseq_results_long = DEseq_results_annotated %>% 
      pivot_longer(cols = UTR3P:CTCF_in_exon, 
                   names_to = "feature", 
                   values_to = "feature_overlap")
    
df_plot = DEseq_results_long %>% 
  filter(feature == "FANTOM5_promoter" & feature_overlap == TRUE | feature == "Intergenic" & feature_overlap == TRUE )
    
plot <- ggplot() +
  geom_point(
    data = df_plot, 
    aes(x = log2FoldChange, y = -1 * log10(pvalue), colour = feature), 
    size = 3, alpha = 0.5, fill = NA, shape = 21, stroke = 2) +
  geom_hline(aes(yintercept = 2), colour = "red", linetype = "dashed") +
  geom_text(label = "p-value = 0.01", colour = "red", aes(x = -6, y = 3)) +
  scale_color_manual(name = "Genomic feature", values = c("#0072B2","#F0E442")) +
  scale_x_continuous(limits = c(-7, 7)) + 
  ylim(NA, 100) +
  labs(subtitle = name,
       x = "log2(FoldChange)",
       y = "-log10(Pvalue)") +
  # theme_tufte()+
  theme(
    axis.line.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    axis.title = element_text(size = 10),
    plot.subtitle = element_text(size = 14, hjust = 0.5))

print(plot)

ggsave(plot = plot, filename = paste0(dir_output, "/volcano_promoter_intergenic_",name,".png"))
    




##############################################################
#### Calcul du fold change ??? :
# cond_avant = trois donneurs = 0 0 0
# cond_après = trois donneurs = 0 7 2
# On fait +1 et ensuite la moyenne des donneurs soit :
# cond_avant = (1+1+1)/3 = 1
# cond_après = (1+8+3)/3 = 4
# Puis on divise après par avant = log2foldchange = 4/1 = 4
#############################################################











#####################################################################################
################# CODE AUTOMATISE POUR PLUSIEURS ECHANTILLONS ####################### => très lent pour les readcounts, faire plutôt avec snakefile
#####################################################################################


# ## Avec un seul bloc d'enregistrement 
# 
# gr_files_list = list.files(path = "D_Analysis/ATAC_differential_accessibility",
#                            pattern = "_ann.gr.rds",
#                            full.names = TRUE)
# 
# for (i in 1:length(gr_files_list)) {
#   
#   #**************************************
#   #* Chargement des fichiers nécessaires
#   #**************************************
#   
#   # Pour faire l'analyse différentielle entre deux conditions, on a besoin : 
#   # - le Grange de l'union des conditions (donneurs fusionnés avec threshold pour chaque condition)
#   name_gr_union = gr_files_list[i]
#   gr_union = readRDS(name_gr_union)
#   
#   cond1 = str_extract(name_gr_union, pattern = "(?<=_vs_)[:alnum:]{2,5}_[:digit:]{2}h(?=_D)")
#   cond2 = str_extract(name_gr_union, pattern = "(?<=differential_accessibility/)[:alnum:]{2,5}_[:digit:]{2}h(?=_D)")
#   name = paste0(cond2, "_vs_", cond1)
#   
#   # - le bam de chaque condition pour tous les donneurs (avant passage au bloc2, un fichier par donneur)
#   name_bam_cond1 = str_subset(list.files(path = "A_Initial_data/merged_analysis/bam", pattern = cond1, full.names = TRUE), 
#                               pattern = "_D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
#   name_bam_cond2 = str_subset(list.files(path = "A_Initial_data/merged_analysis/bam", pattern = cond2, full.names = TRUE), 
#                               pattern = "_D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
#   bam_files = c(name_bam_cond1, name_bam_cond2)
#   
#   # - le Grange de chaque condition en donneurs fusionnés avec threshold
#   gr_cond1 = readRDS(str_subset(list.files(path = "A_Initial_data/merged_analysis/genomic_ranges", pattern = cond1, full.names = TRUE), 
#                                 pattern = "_threshold_[:digit:]{2,3}_ann.gr.rds"))
#   gr_cond2 = readRDS(str_subset(list.files(path = "A_Initial_data/merged_analysis/genomic_ranges", pattern = cond2, full.names = TRUE), 
#                                 pattern = "_threshold_[:digit:]{2,3}_ann.gr.rds"))
# 
#   #**********************
#   # Create readcount matrix
#   #**********************
#   
#     # Mise en forme du Grange en data.frame pour le readcount
#     df = data.frame(gr_union) %>%
#       tibble::rownames_to_column(var = "number") %>%
#       dplyr::mutate(GeneID = paste0("peak_", number), .keep = "unused") %>%
#       dplyr::rename(Chr = seqnames, Start = start, End = end, Strand = strand) %>%
#       dplyr::select(GeneID, Chr, Start, End, Strand)
#     
#     # Calcul du readcount
#     readCount <- featureCounts(files = bam_files,
#                                annot.ext = df,  
#                                isPairedEnd = TRUE,
#                                nthreads = 1,
#                                countChimericFragments = FALSE,
#                                countMultiMappingReads = TRUE)
#     
#     df_annot = readCount$annotation
#     
#     # Enregistrement des readcount
#     print("Save txt...")
#     matrix_count = readCount$counts
#     write.table(matrix_count, 
#                 file = paste0(dir_output, "test_featurecounts_", name, ".txt"),
#                 sep = "\t", quote = FALSE)
#     
#     print("Save csv...")
#     count_df = tibble::rownames_to_column(data.frame(readCount$counts), "GeneID")
#     write.table(count_df, 
#                 file = paste0(dir_output, "test_readcount_", name, ".csv"),
#                 sep = ";", row.names = FALSE)
#     
#     
#     #**********************
#     # Look on the readcount matrix
#     #**********************
#     
#     # Ajout de l'information peak ou pas peak
#     # Lors du readcount, on perd l'info de la présence d'un peak détecté par peak calling ou non 
#     # On va faire des overlap pour voir si les régions extraites de l'union correspondent à un peak dans chaque condition
#     union_gr = makeGRangesFromDataFrame(df = df_readcount, 
#                                         keep.extra.columns = TRUE,
#                                         start.field = "Start",
#                                         end.field = "End",
#                                         seqnames.field = "Chr",
#                                         ignore.strand = FALSE)
#     
#     tab = data.frame(peak_cond1 = as.logical(countOverlaps(union_gr, gr_cond1)),
#                      peak_cond2 = as.logical(countOverlaps(union_gr, gr_cond2))) %>%
#       tibble::rownames_to_column(var = "region") %>%
#       dplyr::mutate(region = paste0("region_", region))
#     
#     df_readcount = left_join(df_readcount, tab, by = "region")
#     
#     
#     # Pour les volcano, on garde seulement les régions qui correspondent à un peak dans les deux conditions 
#     
#     df_peaks_persistant = df_readcount %>% 
#       dplyr::filter(peak_cond1 == TRUE & peak_cond2 == TRUE) %>%
#       dplyr::select(-peak_cond1, -peak_cond2)
#     
#     df_type_peak = data.frame(differential = name,
#                               persistant_peaks = nrow(df_peaks_persistant),
#                               close_open_peaks = nrow(df_readcount %>% dplyr::filter(peak_cond1 == FALSE & peak_cond2 == TRUE)),
#                               open_close_peaks = nrow(df_readcount %>% dplyr::filter(peak_cond1 == TRUE & peak_cond2 == FALSE)))
#     
#     matrix_peaks_volcano = df_peaks_persistant[,7:ncol(df_peaks_persistant)]
#     
#     
#     #**********************
#     # Deseq on readcount matrix
#     #**********************
#     
#     # DEseq2 parameters (chaque sample a 3 donneurs donc 6 lignes au total)
#     coldata <- data.frame(sample = str_extract(bam_files, pattern = "(?<=bam/)[:alnum:]{2,5}_[:digit:]{2}h_D[:digit:]{1}(?=.bam)")) %>%
#       dplyr::mutate(condition = ifelse(str_detect(sample, pattern = cond1), "before", "after"), 
#                     type = "paired-end")
#     
#     # Run DEseq
#     # Dans coldata, une ligne correspond à une colonne de count_df
#     # Normalisation des données avec DEseq
#     dds <- DESeqDataSetFromMatrix(
#       countData = matrix_peaks_volcano,
#       colData = coldata,
#       design = ~ condition)
#     
#     ## Setting the samples tagged "before" as reference
#     dds$condition <- relevel(dds$condition, ref = "before") 
#     dds <- DESeq(dds)
#     
#     # Affichage de la matrice normalisée
#     cm = data.frame(counts(dds, normalized=TRUE))
#     
#     res <- results(dds)
#     res = as_tibble(res)
#     res = res %>% 
#       mutate(region = paste0(name,"_region_", 1:nrow(res))) %>%
#       relocate(region, .before = baseMean)
#     
#     #**********************
#     # Annotate DESeq2 results
#     #**********************
#     
#     # 1 Annotate the region resulting from the union of peaks at two consecutive times ## Idem script bloc 2 : annotate_grange.R
#     
#     all_annotations = loadRData(paste0(dir, "A_Initial_data/Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_woThisto_FANTOM5_prom_gr.rda"))
#     FANTOM_prom_gr = loadRData(paste0(dir, "A_Initial_data/prom_gene_fantom_gr.rdata"))
#     annotations_types = levels(factor(all_annotations$annotation))
#     
#     # First a matrix is created filled with FALSE and added to the Grange
#     metadata = matrix(FALSE, ncol = length(annotations_types), nrow = length(gr_union))
#     colnames(metadata) = annotations_types
#     mcols(gr_union) = metadata
#     
#     # for each of the annotations types an overlap is calculated and used to assigned the peak as TRUE when overlapping with the annotation
#     for (i in 1:ncol(metadata)){
#       sub_annot = all_annotations[all_annotations$annotation == annotations_types[i]]
#       overlaps = findOverlaps(gr_union, sub_annot)
#       mcols(gr_union)[queryHits(overlaps),i] = TRUE
#     }
#     
#     # Adding overlap with fantom5 promoter database
#     overlap_with_FANTOM = findOverlaps(gr_union, FANTOM_prom_gr)
#     mcols(gr_union)[,"FANTOM5_promoter"] = FALSE
#     mcols(gr_union)[queryHits(overlap_with_FANTOM), "FANTOM5_promoter"] = TRUE
#     
#     colnames(mcols(gr_union)) = c("UTR3P","UTR5P","CpG", "CTCF","Exons","Introns","TSS_mp1kb","FANTOM5_promoter")
#     
#     mcols(gr_union) = as_tibble(mcols(gr_union)) %>%
#       dplyr::mutate(Intergenic = ifelse(UTR3P == FALSE & UTR5P == FALSE & Exons == FALSE & Introns == FALSE & FANTOM5_promoter == FALSE & TSS_mp1kb == FALSE, TRUE, FALSE)) %>%
#       dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
#       dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
#       dplyr::mutate(CTCF_Intergenic = ifelse(Intergenic == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
#       dplyr::mutate(CTCF_in_intron = ifelse(Introns == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
#       dplyr::mutate(CTCF_in_exon = ifelse(Exons == TRUE & CTCF == TRUE, TRUE, FALSE))
#     
#     df = as_tibble(gr_union) %>% 
#       tibble::rownames_to_column(var = "number") %>%
#       dplyr::mutate(region = paste0(name, "_region_", number)) %>%
#       dplyr::select(-number) %>%
#       relocate(region, .before = seqnames)
#     
#     # 2 Adding the information concerning region resulting from the union to DEseq results
#     
#     DEseq_results = res %>% 
#       mutate(regulation = case_when(
#         pvalue < 0.01 ~ "significative",
#         pvalue > 0.01 ~ "Non-significative"))
#     
#     # Adding start|end|seqnames informations to regions
#     DEseq_results_annotated = left_join(DEseq_results,  df, by = "region") 
#     
#     ## enlever les NA ! 
#     
#     write.csv2(x = DEseq_results_annotated, file = paste0(dir_output, "/DEseq_results_", name, ".csv"))
#     
#     df_recap_peak = table(DEseq_results_annotated$regulation)   
#     
#     
#     #**********************
#     # Volcano_plot
#     #**********************
#     
#     DEseq_results_long = DEseq_results_annotated %>% 
#       pivot_longer(cols = UTR3P:CTCF_in_exon, 
#                    names_to = "feature", 
#                    values_to = "feature_overlap")
#     
#     df_plot = DEseq_results_long %>% 
#       filter(feature == "FANTOM5_promoter" & feature_overlap == TRUE | feature == "Intergenic" & feature_overlap == TRUE )
#     
#     plot <- ggplot() +
#       geom_point(
#         data = df_plot, 
#         aes(x = log2FoldChange, y = -1 * log10(pvalue), colour = feature), 
#         size = 3, alpha = 0.5, fill = NA, shape = 21, stroke = 2) +
#       geom_hline(aes(yintercept = 2), colour = "red", linetype = "dashed") +
#       geom_text(label = "p-value = 0.01", colour = "red", aes(x = -6, y = 3)) +
#       scale_color_manual(name = "Genomic feature", values = c("#0072B2","#F0E442")) +
#       scale_x_continuous(limits = c(-7, 7)) + 
#       ylim(NA, 100) +
#       labs(subtitle = name,
#            x = "log2(FoldChange)",
#            y = "-log10(Pvalue)") +
#       # theme_tufte()+
#       theme(
#         axis.line.y = element_line(color = "black"),
#         axis.line.x = element_line(color = "black"),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 11),
#         axis.title = element_text(size = 10),
#         plot.subtitle = element_text(size = 14, hjust = 0.5))
#     
#     print(plot)
#     
#     ggsave(plot = plot, filename = paste0(dir_output, "/volcano_promoter_intergenic_",name,".png"))
#     
# }
# 
# 
