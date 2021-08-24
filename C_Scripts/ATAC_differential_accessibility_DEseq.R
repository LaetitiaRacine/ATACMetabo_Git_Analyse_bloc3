
#**********************
# Command line to call the script in linux consol
#**********************

" DEseq : étude différentielle des peaks présents dans les deux conditions testées (Les peaks sont-ils plus grands ou plus petits ?)

Usage:
  ATAC_differential_accessibility.R [options] create_readcountmatrix <bam_dir> <gr_dir> <gr_union> <readcount_csv>
  ATAC_differential_accessibility.R [options] deseq_object <readcount_csv> <gr_annot_union> <norm_matrix> <deseq_result>
  ATAC_differential_accessibility.R [options] annot_volcano <deseq_result> <plot_name>
  ATAC_differential_accessibility.R -h | --help

Options:
  --output_txt <matrix_txt>           Matrix count as txt file
  --output_peakstype <peakstype_csv>  Dataframe of peaks categories as csv file
  -h, --help                          Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#**********************
# Libraries and function loading
#**********************

suppressPackageStartupMessages({
  suppressWarnings(library(GenomicRanges))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(DESeq2))
  suppressWarnings(library(stringr))
  suppressWarnings(library(tibble))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(Rsubread))
})

#*********************************************************************************************************
#* Pour créer la matrice readcount et ajouter les informations sur les peaks 
#* *******************************************************************************************************

if (arguments$create_readcountmatrix) { 
  
  #********************
  ### Inputs' loading
  #********************
 
  # Grange de l'union des conditions (donneurs fusionnés avec threshold pour chaque condition)
    name_gr_union = arguments$gr_union
    gr_union = readRDS(name_gr_union)
  
    cond1 = str_extract(name_gr_union, pattern = "(?<=_vs_)[:alnum:]{2,5}_[:digit:]{2}h(?=_D)")
    cond2 = str_extract(name_gr_union, pattern = "(?<=differential_accessibility/)[:alnum:]{2,5}_[:digit:]{2}h(?=_D)")
    name = paste0(cond2, "_vs_", cond1)
  
    # bam de chaque condition pour tous les donneurs (avant passage au bloc2, un fichier par donneur)
    name_bam_cond1 = str_subset(list.files(path = arguments$bam_dir, pattern = cond1, full.names = TRUE), pattern = "_D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
    name_bam_cond2 = str_subset(list.files(path = arguments$bam_dir, pattern = cond2, full.names = TRUE), pattern = "_D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
    bam_files = c(name_bam_cond1, name_bam_cond2)
    
    # Grange de chaque condition en donneurs fusionnés avec threshold
    gr_cond1 = readRDS(str_subset(list.files(path = arguments$gr_dir, pattern = cond1, full.names = TRUE), pattern = "_threshold_[:digit:]{2,3}_ann.gr.rds"))
    gr_cond2 = readRDS(str_subset(list.files(path = arguments$gr_dir, pattern = cond2, full.names = TRUE), pattern = "_threshold_[:digit:]{2,3}_ann.gr.rds"))
  
  #**********************
  # Create readcount matrix
  #**********************
  
    # Mise en forme du Grange en data.frame pour le readcount
    df_gr_union = data.frame(gr_union) %>%
      tibble::rownames_to_column(var = "number") %>%
      dplyr::mutate(GeneID = paste0("peak_", number), .keep = "unused") %>%
      dplyr::rename(Chr = seqnames, Start = start, End = end, Strand = strand) %>%
      dplyr::select(GeneID, Chr, Start, End, Strand)
  
    # Calcul du readcount
    readCount <- featureCounts(files = bam_files,
                               annot.ext = df_gr_union,  
                               isPairedEnd = TRUE,
                               nthreads = 1,
                               countChimericFragments = FALSE,
                               countMultiMappingReads = TRUE)
  
    # Ajout de l'information peak ou pas peak dans le peak calling initial
    df_readcount = left_join(readCount$annotation, tibble::rownames_to_column(data.frame(readCount$counts), "GeneID"), by = "GeneID") %>%
      dplyr::mutate(region = str_replace(GeneID, "peak", "region"), .before = Chr) %>%
      dplyr::select(-GeneID) %>%
      dplyr::mutate(samples = name, .before = region )

    gr_readcount = makeGRangesFromDataFrame(df = df_readcount, keep.extra.columns = TRUE, ignore.strand = FALSE,
                                            start.field = "Start", end.field = "End", seqnames.field = "Chr")
    
    correspondance = data.frame(peak_cond1 = as.logical(countOverlaps(gr_readcount, gr_cond1)),
                                peak_cond2 = as.logical(countOverlaps(gr_readcount, gr_cond2))) %>%
      tibble::rownames_to_column(var = "region") %>%
      dplyr::mutate(region = paste0("region_", region))
                                        
    df_readcount = left_join(df_readcount, correspondance, by = "region")
    
    # Classement des peaks dans des catégories
    df_type_peak = data.frame(differential = name,
                              persistant_peaks = nrow(df_readcount %>% dplyr::filter(peak_cond1 == TRUE & peak_cond2 == TRUE)),
                              close_open_peaks = nrow(df_readcount %>% dplyr::filter(peak_cond1 == FALSE & peak_cond2 == TRUE)),
                              open_close_peaks = nrow(df_readcount %>% dplyr::filter(peak_cond1 == TRUE & peak_cond2 == FALSE)))
    
    #**********************
    # Outputs' saving
    #**********************
    
    print("Save readcount file with peak information in csv format.")
    write.table(df_readcount, file = arguments$readcount_csv, sep = ";", row.names = FALSE)
    
    if (!is.null(arguments$output_txt)) {
      print("Save matrix count file in txt format.")
      matrix_count = readCount$counts 
      rownames(matrix_count) <- str_replace(rownames(matrix_count),"peak", "region")
      write.table(matrix_count, file = arguments$output_txt, sep = "\t", quote = FALSE)
    }
    
    if (!is.null(arguments$output_peakstype)) {
      print("Save table with peaks categories in csv format.")
      write.table(df_type_peak, file = arguments$output_peakstype, sep = ";", row.names = FALSE )
    }
    
}

#*********************************************************************************************************
#* Pour normaliser les données et calucler les statistiques DEseq
#* *******************************************************************************************************

if (arguments$deseq_object) {
  
  #********************
  ### Inputs' loading
  #********************

  df_readcount = read.csv2(arguments$readcount_csv)
  cond1 = str_extract(string = unique(df_readcount$samples), pattern = "(?<=_vs_).+")
  cond2 = str_extract(string = unique(df_readcount$samples), pattern = ".+(?=_vs_)")
  name = paste0(cond2, "_vs_", cond1)
  
  gr_union = readRDS(arguments$gr_annot_union)
  
  #**********************
  # Deseq on readcount matrix
  #**********************
  
  # Pour le deseq et les volcano, on garde seulement les régions qui correspondent à un peak dans les deux conditions 
  matrix_peaks = df_readcount %>% 
    dplyr::filter(peak_cond1 == TRUE & peak_cond2 == TRUE) %>%
    dplyr::select(-peak_cond1, -peak_cond2, -samples, -Chr, -Start, - End, -Strand, -Length) %>%
    column_to_rownames(var = "region")
  
  # DEseq2 parameters => experience design data frame (one row of coldata correspond to one column of matrix_peaks)
  coldata <- data.frame(sample = colnames(matrix_peaks)) %>%
    dplyr::mutate(condition = ifelse(str_detect(sample, pattern = cond1), "before", "after")) %>%
    dplyr::mutate(type = "paired-end")
  
  # Normalize data with DEseq
  # Dans coldata, une ligne correspond à une colonne de count_df
  dds <- DESeqDataSetFromMatrix(
    countData = matrix_peaks,
    colData = coldata,
    design = ~ condition)
  
  ## Setting the samples tagged "before" as reference
  dds$condition <- relevel(dds$condition, ref = "before") 
  dds <- DESeq(dds)
  
  # Normalized matrix
  cm = data.frame(counts(dds, normalized=TRUE))
  
  # Statistics from deseq
  res = as.tibble(results(dds)) 
  res = res %>% 
    dplyr::mutate(region = rownames(cm), .before = baseMean) %>%
    dplyr::mutate(regulation = case_when(pvalue < 0.01 ~ "significative", pvalue > 0.01 ~ "Non-significative")) 

  # Add annotation in deseq object
  df_union = as_tibble(gr_union) %>%
    tibble::rownames_to_column(var = "number") %>%
    dplyr::mutate(sample = name, .before = seqnames) %>%
    dplyr::mutate(region = paste0("region_", number), .before = seqnames) %>%
    dplyr::select(-number) 
  
  # Keep only regions present in DEseq_results and add information from df_union
  DEseq_results_annotated = right_join(df_union, res, by = "region") 
  print(table(DEseq_results_annotated$regulation))   
  
  #**********************
  # Outputs' saving
  #**********************
  
  print("Save normalized matrix count in txt format.")
  write.table(cm, file = arguments$norm_matrix, sep = "\t",  quote = FALSE)
  
  print("Save deseq statistics in csv format.")
  write.table(DEseq_results_annotated, file = arguments$deseq_result, sep = ";", row.names = FALSE)

}

#*********************************************************************************************************
#* Pour tracer les volcano plots
#* *******************************************************************************************************

if (arguments$annot_volcano) {

  #********************
  ### Inputs' loading
  #********************

  DEseq_results_annotated = read.csv2(file = arguments$deseq_result, sep =";")
  start = grep("baseMean", colnames(DEseq_results_annotated))
  end = grep("padj", colnames(DEseq_results_annotated))
  DEseq_results_annotated[,start:end] = lapply(DEseq_results_annotated[,start:end], as.numeric)
  name = unique(DEseq_results_annotated$sample)
  
  #********************
  ### Create volcano plots
  #********************

  # Pivot du tableau pour permettre de tracer le plot  
  DEseq_results_long = DEseq_results_annotated %>% 
    pivot_longer(cols = UTR3P:CTCF_in_exon, 
                 names_to = "feature", 
                 values_to = "feature_overlap") 
  
  # On conserve seulement les peaks qui sont dans les promoters et dans l'intergenique
  df_plot = DEseq_results_long %>% 
    filter(feature == "FANTOM5_promoter" & feature_overlap == TRUE | feature == "Intergenic" & feature_overlap == TRUE )
  
  plot <- ggplot() +
    geom_point(data = df_plot, aes(x = log2FoldChange, y = -1 * log10(pvalue), colour = feature), size = 3, alpha = 0.5, fill = NA, shape = 21, stroke = 2) +
    geom_hline(aes(yintercept = 2), colour = "red", linetype = "dashed") +
    geom_text(label = "p-value = 0.01", colour = "red", aes(x = -6, y = 3)) +
    scale_color_manual(name = "Genomic feature", values = c("#0072B2","#F0E442")) +
    scale_x_continuous(limits = c(-7, 7)) + 
    ylim(NA, 100) +
    labs(subtitle = name, x = "log2(FoldChange)", y = "-log10(Pvalue)") +
    # theme_tufte()+
    theme(axis.line.y = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      axis.title = element_text(size = 10),
      plot.subtitle = element_text(size = 14, hjust = 0.5))

  ggsave(plot = plot, 
         filename = arguments$plot_name)

}
