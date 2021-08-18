
dir = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/"
dir_output = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/D_Analysis/ATAC_differential_accessibility"

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls()!="fileName"])
}

### DESeq2 on readcount matrix (interval union)

matrix_readcount_union_files = list.files(path = dir_output,
                                          pattern = "readcount.+.csv",
                                          full.names = TRUE)
DEseq_res_list = list()

for (i in 1:length(matrix_readcount_union_files)) {
  
  matrix_count_union = read.csv2(matrix_readcount_union_files[i], sep = ";")
  name_vs = str_extract(matrix_readcount_union_files[i], pattern = "(?<=readcount_).+(?=.csv)")
  
  # DEseq2 parameters
  coldata <- data.frame(
    condition = c("before", "after"), 
    type = rep("paired-end", 2)) 
  
  # Run DEseq
  dds <- DESeqDataSetFromMatrix(
    countData = matrix_count_union[,-1], 
    colData = coldata, 
    design = ~ 1)
  
  ## Setting the samples tagged "before" as reference
  # dds$condition <- relevel(dds$condition, ref = "before") 
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  res = as_tibble(res)
  res = res %>% 
    mutate(region = paste0(name_vs,"_region_", 1:nrow(res))) %>%
    relocate(region, .before = baseMean)
  
  DEseq_res_list[[i]] = res
  
}


### Annotate DESeq2 results


# 1 Annotate the region resulting from the union of peaks at two consecutive times ## Idem script bloc 2 : annotate_grange.R

all_annotations = loadRData(paste0(dir, "A_Initial_data/Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_woThisto_FANTOM5_prom_gr.rda"))
FANTOM_prom_gr = loadRData(paste0(dir, "A_Initial_data/prom_gene_fantom_gr.rdata"))
annotations_types = levels(factor(all_annotations$annotation))

peaks_union = list.files(path = dir_output, pattern = "gr.rds", full.names = TRUE)
union_region_list = list()
for (i in 1:length(peaks_union)) {
  name = str_extract(peaks_union[[i]], pattern = "(?<=accessibility/).+(?=_threshold)")
  union_region_list[[i]] = readRDS(peaks_union[[i]])
  names(union_region_list)[i] = name
}


for(interval in 1:length(union_region_list)){
  
  gr = union_region_list[[interval]]
  name = names(union_region_list)[interval]
  
  # First a matrix is created filled with FALSE and added to the Grange
  metadata = matrix(FALSE, ncol = length(annotations_types), nrow = length(gr))
  colnames(metadata) = annotations_types
  mcols(gr) = metadata

  # for each of the annotations types an overlap is calculated and used to assigned the peak as TRUE when overlapping with the annotation
  for (i in 1:ncol(metadata)){
    sub_annot = all_annotations[all_annotations$annotation == annotations_types[i]]
    overlaps = findOverlaps(gr, sub_annot)
    mcols(gr)[queryHits(overlaps),i] = TRUE
  }

  # Adding overlap with fantom5 promoter database
  overlap_with_FANTOM = findOverlaps(gr, FANTOM_prom_gr)
  mcols(gr)[,"FANTOM5_promoter"] = FALSE
  mcols(gr)[queryHits(overlap_with_FANTOM), "FANTOM5_promoter"] = TRUE

  colnames(mcols(gr)) = c("UTR3P","UTR5P","CpG", "CTCF","Exons","Introns","TSS_mp1kb","FANTOM5_promoter")

  mcols(gr) = as_tibble(mcols(gr)) %>%
    dplyr::mutate(Intergenic = ifelse(UTR3P == FALSE & UTR5P == FALSE & Exons == FALSE & Introns == FALSE & FANTOM5_promoter == FALSE & TSS_mp1kb == FALSE, TRUE, FALSE)) %>%
    dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
    dplyr::mutate(CpG_Intergenic = ifelse(Intergenic == TRUE & CpG == TRUE, TRUE, FALSE)) %>%
    dplyr::mutate(CTCF_Intergenic = ifelse(Intergenic == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
    dplyr::mutate(CTCF_in_intron = ifelse(Introns == TRUE & CTCF == TRUE, TRUE, FALSE)) %>%
    dplyr::mutate(CTCF_in_exon = ifelse(Exons == TRUE & CTCF == TRUE, TRUE, FALSE))

  df = as_tibble(gr) %>% 
    tibble::rownames_to_column(var = "number") %>%
    dplyr::mutate(region = paste0(name, "_region_", number)) %>%
    dplyr::select(-number) %>%
    relocate(region, .before = seqnames)
  
  union_region_list[[interval]] = df
  names(union_region_list)[interval] = paste0(name,"_gr_annotated")

}  
  
# 2 Adding the information concerning region resulting from the union to DEseq results

for(i in 1:(length(DEseq_res_list))){
  
  DEseq_results = DEseq_res_list[[i]]
  
  DEseq_results = DEseq_results %>% 
    mutate(
      regulation = case_when(
        pvalue < 0.01 ~ "significative",
        pvalue > 0.01 ~ "Non-significative"))
  
  # Adding start|end|seqnames informations to regions
  colnames(union_region_list[[i]])[1] = "region"
  
  DEseq_results_annotated = left_join(
    DEseq_results,
    union_region_list[[i]],
    by = "region") 
  
  write.csv2(x = DEseq_results_annotated, file = paste0(dir_output, "/DEseq_results_", names(union_region_list)[[i]], ".csv"))
}
  
  
### Volcano plot

DEseq_files = list.files(path = dir_output, pattern = "DEseq_results", full.names = T)

for(i in 1:length(DEseq_files)){
  
  DEseq_results = read.csv2(DEseq_files[i])[,-1] # pour enlever la colonne X qui se rajoutait toute seule
  
  DEseq_results_long = DEseq_results %>% 
    pivot_longer(
      cols = UTR3P:CTCF_in_exon, 
      names_to = "feature", 
      values_to = "feature_overlap")
  
  df_plot = DEseq_results_long %>% 
    filter(feature == "FANTOM5_promoter" & feature_overlap == TRUE | feature == "Intergenic" & feature_overlap == TRUE )
  
  name = str_extract(df_plot$region[1], pattern = ".+(?=_region)")
  
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
  
  # print(plot)

  ggsave(plot = plot, filename = paste0(dir_output, "/volcano_promoter_intergenic_",name,".png"))
  
} 

















  
  
# #############
# ### 3 : Differential analysis between control (t) and treated (t+1) datasets
# ###     (DEseq2, foldChange and pValue associated calculation)
# #############
# 
# # Only keep number of reads detected for each samples in the union region
# matrix_count = readCount$counts
# #colnames(matrix_count) <- unlist(conditions_list[cond_id:(cond_id+1)])
# colnames(matrix_count) <- samples
# 
# write.table(matrix_count, file = paste0(the_dir, "results/readCount_output/featuresCount_",
#                                         name_condition_list_b, "_vs_", name_condition_list_a,
#                                         ".txt"), sep = "\t", quote = FALSE)
# 
# print("Fonction read_deseq partie 3 - matrice ok")
# print(paste0("Output = featuresCount_", name_condition_list_b, "_vs_", name_condition_list_a,".txt"))
# 
# ### Parameters  for DEseq
# # seqType <- rep("paired-end",6)
# seqType <- rep("paired-end", length(samples))
# # condition <- c(rep("control",3),rep("treated",3)) #Control = t | Treated = t+1
# condition <- c(rep("control", (length(samples)/2)), rep("treated", (length(samples)/2)))
# # on part du principe que pour comparer deux time_point, on doit avoir le même nombre de sample dans chaque
# # coldata <- as.data.frame(cbind(condition = condition,  type = seqType), row.names = unlist(conditions_list[cond_id:(cond_id+1)]) )
# coldata <- as.data.frame(cbind(condition = condition,  type = seqType), row.names = samples)
# 
# # Prepare DEseq object
# dds <- DESeqDataSetFromMatrix(matrix_count, colData = coldata, design = ~ condition)
# 
# # Filter peaks base on minimum reads
# keep <- rowSums(counts(dds) >= 10) >= 2 # At least 2 samples must have a minimum of 10 reads to consider the region
# dds <- dds[keep,]
# 
# # Calculate DEseq Matrix
# dds <- DESeq(dds)
# res <- results(dds)
# 
# # Merging peaks informations and DEseq results
# colnames(differential_peaks_union_df)[1] <- "name"
# d2res <- as.data.frame(res)
# d2res$name <- rownames(d2res)
# d2res <- merge(differential_peaks_union_df[,-5], d2res, by=c("name")) # Adding chr/start/end informations
# if (sum(is.na(d2res$pvalue)) != 0) d2res = d2res[-which(is.na(d2res$pvalue)),]   ### Remove rows where DEseq generated NA for pValues
# 
# ## Indicate differential status of the region between 2 conditions
# d2res = d2res %>%
#   mutate(regulation = case_when(log2FoldChange < 0 ~ "gained-close",
#                                 log2FoldChange > 0 ~ "gained-open"))
# # Export final results
# write.table(d2res, file = paste0(the_dir, "results/DEseq_output/DEseq2_results_",
#                                  name_condition_list_b, "_vs_", name_condition_list_a,
#                                  ".txt"), sep = "\t", quote = FALSE)
# 
# print ("Fonction read_deseq partie 3 OK")
# print(paste0("Output = DEseq2_result_", name_condition_list_b, "_vs_", name_condition_list_a, ".txt"))
# 
# ############
# ### 4: Create genomic.ranges from DEseq2 results
# ############
# 
# Peak_Grange = makeGRangesFromDataFrame(d2res[c(1,2,3,4,11)], # Keeping only peak info + regulation
#                                        keep.extra.columns = TRUE,
#                                        ignore.strand = FALSE,
#                                        seqinfo = NULL,
#                                        seqnames.field= c ("seqnames", "seqname",
#                                                           "chromosome", "chrom",
#                                                           "chr", "chromosome_name",
#                                                           "seqid"),
#                                        start.field = "start",
#                                        end.field = c("end", "stop"),
#                                        strand.field = "strand",
#                                        starts.in.df.are.0based = FALSE)
# 
# # Create and save increasing regions
# peaks_inc_gr = Peak_Grange[Peak_Grange$regulation == "gained-open"]
# peaks_inc_gr = clean_gr(peaks_inc_gr)
# 
# print("Fonction read_deseq partie 4 increasing ok")
# print(paste0("Output = peaks_inc_", name_condition_list_b, "_vs_", name_condition_list_a, "_gr.rda"))
# 
# save(peaks_inc_gr,
#      file = paste0(the_dir,"results/genomic_ranges/differential_peaks/peaks_inc_",
#                    name_condition_list_b, "_vs_", name_condition_list_a, "_gr.rda"))
# 
# 
# # Create and save decreasing regions
# peaks_dec_gr = Peak_Grange[Peak_Grange$regulation == "gained-close"]
# peaks_dec_gr = clean_gr(peaks_dec_gr)
# 
# print("Fonction read_deseq partie 4 decreasing ok")
# print(paste0("Output = peaks_dec_", name_condition_list_b, "_vs_", name_condition_list_a, "_gr.rda"))
# 
# save(peaks_dec_gr,
#      file = paste0(the_dir,"results/genomic_ranges/differential_peaks/peaks_dec_",
#                    name_condition_list_b, "_vs_", name_condition_list_a, "_gr.rda"))
# }
# 
# 
