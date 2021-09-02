
# Gene coverage with Gviz
# https://www.bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html 

#########################
### Command line to call the script in linux consol
#########################

"Create plot for gene coverage

Usage:
  ATAC_gene_coverage.R [options] multi <prom_file> <region_file> <bam_dir> <gr_dir> 
  ATAC_gene_coverage.R [options] single <prom_file> <region_file> <bam_dir> <peak_file>
  ATAC_gene_coverage.R [options] multi_choice <prom_file> <region_file> <bam_dir> <peak_file>...
  ATAC_gene_coverage.R -h | --help

Options:
  -h, --help              Show this screen

" -> doc

library(docopt)
arguments <- docopt(doc)

#########
### Initialization : libraries loading and function definitions
#########

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(GenomicRanges))
  suppressWarnings(library(Gviz))
  suppressWarnings(library(stringr))
  
})

col_fun = function(file_name) {
  color = case_when(str_detect(file_name, "2DG") ~ "#F8766D",
                    str_detect(file_name, "DON") ~ "#00C094",
                    str_detect(file_name, "aK") ~ "#C49A00",
                    str_detect(file_name, "VPA") ~ "#A58AFF",
                    str_detect(file_name, "AOA") ~"#53B400",
                    str_detect(file_name, "MP") ~ "#00B6EB",
                    str_detect(file_name, "Xvivo") ~ "#FB61D7")
  return(color)
}


####################
## Loading of files
####################

prom_gene_fantom_gr = readRDS(arguments$prom_file)
df_region = read.csv2(file = arguments$region_file, sep = ";")

axTrack <- GenomeAxisTrack() # scale
blank_space = GenomeAxisTrack(col = "white",fill = "white", fontcolor = "white")   # to improve plot readibility

##############################
# One plot with all conditions 
##############################

if (arguments$multi == TRUE) {
  
#########
# Loading and storing files
#########

print("Chargement des fichiers")
  
# Extract peaks files and bam directions
list_peaks_files = list.files(arguments$gr_dir, pattern = "threshold_10_ann.gr.rds", full.names = TRUE)
names(list_peaks_files) = str_extract(string = list_peaks_files, pattern = "(?<=/genomic_ranges/).+(?=_threshold_10_ann.gr.rds)")
list_bam_files = str_subset(list.files(arguments$bam_dir, pattern = ".bam", full.names = TRUE), pattern = "-D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
names(list_bam_files) = str_extract(string = list_bam_files, pattern = "(?<=bam/).+(?=.bam)")

# Load peak files in a list and name it
list_peaks_kept = lapply(X = list_peaks_files, FUN = readRDS)

#########
# Tracking plots creation for each gene
#########

for (i in 1:nrow(df_region)) { # one plot per gene of the list
   
  print(paste("Boucle : gène", df_region$gene[i]))
  
  gr_region = GRanges(seqnames = df_region$chr[i], 
                      ranges= IRanges(start = df_region$start[i], 
                                      end = df_region$end[i]))
  
  # Chromosom with gene illustration
  # ********************************
  axTrack <- GenomeAxisTrack() # scale
  idxTrack <- IdeogramTrack(genome="hg19", chromosome= df_region$chr[i]) # drawing
  blank_space = GenomeAxisTrack(col = "white",fill = "white", fontcolor = "white")   # to improve plot readibility
  
  # Isoform tracking for the studied gene and promoters
  # ***************************************************
  ensGenes <- UcscTrack(genome = "hg19", 
                        chromosome = df_region$chr[i],
                        background.title = "transparent",
                        alpha.title = 1, 
                        col.title = "black",
                        cex.title = 1.5, 
                        track = "ensGene",
                        from = df_region$start[i], 
                        to = df_region$end[i],
                        trackType = "GeneRegionTrack", 
                        rstarts = "exonStarts",
                        rends = "exonEnds", 
                        transcript = "name", 
                        strand = "strand",
                        fill = "lightgrey", 
                        name = df_region$gene[i])
  
  # Promoteurs en vert sur le track des isoformes
  # **********************************************

  prom_peak = subsetByOverlaps(prom_gene_fantom_gr, gr_region)
  
  df_prom_retained = tibble("start" = start(prom_peak),
                            "end" = end(prom_peak),
                            "chr" = as.character(seqnames(prom_peak)))
  
  if(nrow(df_prom_retained) != 0){ 
    promoterTrack = HighlightTrack(trackList = ensGenes,
                                   start = df_prom_retained$start,
                                   end = df_prom_retained$end,
                                   chromosome = df_region$chr[i],
                                   col = "green", 
                                   fill = "green")
  }
  
  
  # Regroupement des pistes dans une liste de plots à tracer 
  # ********************************************************
  
  if (exists("promoterTrack") == TRUE){
    list_to_plot = list(idxTrack, axTrack, blank_space, promoterTrack)
  } else{
    list_to_plot = list(idxTrack, axTrack, blank_space, ensGenes)
  }
  
  print("Datatracking")
        
  # DataTracking 
  # *************
  
  list_datatrack = list()
  
  for (j in 1:length(list_peaks_kept)) { 
    
    name = names(list_peaks_kept)[j]
    
    dt = DataTrack(type = "histogram",
                   name = names(list_peaks_kept[j]),
                   background.title = "transparent",
                   fill.histogram = "#0072B2",
                   col.histogram = "#0072B2",
                   alpha.title = 1,
                   col.title = col_fun(name),
                   col.axis = "black",
                   cex.title = 1,
                   range = unlist(list_bam_files[name]),
                   genome = "hg19",
                   ylim = c(0, df_region$max_read[i]),
                   window = -1,
                   chromosome = df_region$chr[i]
    )
    
    list_datatrack[j] = dt
    names(list_datatrack)[j] = name
    
  }
  
  
  for (k in 1:length(list_peaks_kept)) { 
    
    name = names(list_peaks_kept)[k]
    
    # Determine if peaks are detected in this gene for this time_point
    peaks_in = subsetByOverlaps(list_peaks_kept[[name]],
                                GRanges(seqnames = df_region$chr[i],
                                        ranges = IRanges(start = df_region$start[i], end = df_region$end[i])))
    
    if (length(peaks_in) == 0){  # if not detected, no graph for this time_point for this gene
      
      print(paste0(name, "no_peaks_detected"))
      
    } 
    else{  # if detected, graph with highligted peaks
      
      df_peaks_retained = tibble("start" = start(peaks_in),
                                 "end" = end(peaks_in),
                                 "chr" = as.character(seqnames(peaks_in)))
      
      ht = HighlightTrack(trackList = list_datatrack[name],
                          start = df_peaks_retained$start,
                          end = df_peaks_retained$end,
                          chromosome = df_region$chr[i])
      
      list_to_plot = c(list_to_plot, blank_space, ht) 
      
    }
    
  }
  
  # Tracé des plots
  # ***************
  
  pdf(file =  paste0("D_Analysis/ATAC_gene_coverage/", df_region$gene[i], "_Coverage_plot.pdf"), width = 20, height = 40)
  plotTracks(list_to_plot, from = df_region$start[i], to = df_region$end[i], showTitle = TRUE)
  dev.off()
  
}

}

##############################
# One plot with multiple choosen conditions 
##############################

if (arguments$multi_choice == TRUE) {
  
  #########
  # Loading and storing files
  #########
  
  # Extract peaks files and bam directions
  list_peaks_files = arguments$peak_file
  names(list_peaks_files) = str_extract(string = list_peaks_files, pattern = "(?<=/genomic_ranges/).+(?=_threshold_10_ann.gr.rds)")
  list_bam_files = str_subset(list.files(arguments$bam_dir, pattern = ".bam", full.names = TRUE), pattern = "-D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
  names(list_bam_files) = str_extract(string = list_bam_files, pattern = "(?<=bam/).+(?=.bam)")
  
  # Load peak files in a list and name it
  list_peaks_kept = lapply(X = list_peaks_files, FUN = readRDS)
  
  #########
  # Tracking plots creation for each gene
  #########
  
  for (i in 1:nrow(df_region)) { # one plot per gene of the list
    
    gr_region = GRanges(seqnames = df_region$chr[i], 
                        ranges= IRanges(start = df_region$start[i], 
                                        end = df_region$end[i]))
    
    # Chromosom with gene illustration
    # ********************************
    axTrack <- GenomeAxisTrack() # scale
    idxTrack <- IdeogramTrack(genome="hg19", chromosome= df_region$chr[i]) # drawing
    blank_space = GenomeAxisTrack(col = "white",fill = "white", fontcolor = "white")   # to improve plot readibility
    
    # Isoform tracking for the studied gene and promoters
    # ***************************************************
    ensGenes <- UcscTrack(genome = "hg19", 
                          chromosome = df_region$chr[i],
                          background.title = "transparent",
                          alpha.title = 1, 
                          col.title = "black",
                          cex.title = 1.5, 
                          track = "ensGene",
                          from = df_region$start[i], 
                          to = df_region$end[i],
                          trackType = "GeneRegionTrack", 
                          rstarts = "exonStarts",
                          rends = "exonEnds", 
                          transcript = "name", 
                          strand = "strand",
                          fill = "lightgrey", 
                          name = df_region$gene[i])
    
    # Promoteurs en vert sur le track des isoformes
    # **********************************************
    
    prom_peak = subsetByOverlaps(prom_gene_fantom_gr, gr_region)
    
    df_prom_retained = tibble("start" = start(prom_peak),
                              "end" = end(prom_peak),
                              "chr" = as.character(seqnames(prom_peak)))
    
    if(nrow(df_prom_retained) != 0){ 
      promoterTrack = HighlightTrack(trackList = ensGenes,
                                     start = df_prom_retained$start,
                                     end = df_prom_retained$end,
                                     chromosome = df_region$chr[i],
                                     col = "green", 
                                     fill = "green")
    }
    
    
    # Regroupement des pistes dans une liste de plots à tracer 
    # ********************************************************
    
    if (exists("promoterTrack") == TRUE){
      list_to_plot = list(idxTrack, axTrack, blank_space, promoterTrack)
    } else{
      list_to_plot = list(idxTrack, axTrack, blank_space, ensGenes)
    }
    
    
    # DataTracking 
    # *************
    
    list_datatrack = list()
    
    for (j in 1:length(list_peaks_kept)) { 
      
      name = names(list_peaks_kept)[j]
      
      dt = DataTrack(type = "histogram",
                     name = names(list_peaks_kept[j]),
                     background.title = "transparent",
                     fill.histogram = "#0072B2",
                     col.histogram = "#0072B2",
                     alpha.title = 1,
                     col.title = col_fun(name),
                     col.axis = "black",
                     cex.title = 1,
                     range = unlist(list_bam_files[name]),
                     genome = "hg19",
                     ylim = c(0, df_region$max_read[i]),
                     window = -1,
                     chromosome = df_region$chr[i]
      )
      
      list_datatrack[j] = dt
      names(list_datatrack)[j] = name
      
    }
    
    
    for (k in 1:length(list_peaks_kept)) { 
      
      name = names(list_peaks_kept)[k]
      
      # Determine if peaks are detected in this gene for this time_point
      peaks_in = subsetByOverlaps(list_peaks_kept[[name]],
                                  GRanges(seqnames = df_region$chr[i],
                                          ranges = IRanges(start = df_region$start[i], end = df_region$end[i])))
      
      if (length(peaks_in) == 0){  # if not detected, no graph for this time_point for this gene
        
        print(paste0(name, "no_peaks_detected"))
        
      } 
      else{  # if detected, graph with highligted peaks
        
        df_peaks_retained = tibble("start" = start(peaks_in),
                                   "end" = end(peaks_in),
                                   "chr" = as.character(seqnames(peaks_in)))
        
        ht = HighlightTrack(trackList = list_datatrack[name],
                            start = df_peaks_retained$start,
                            end = df_peaks_retained$end,
                            chromosome = df_region$chr[i])
        
        list_to_plot = c(list_to_plot, blank_space, ht) 
        
      }
      
    }
    
    # Tracé des plots
    # ***************
    name = paste(names(list_peaks_files) , collapse = '_')
    pdf(file =  paste0("D_Analysis/ATAC_gene_coverage/", name, "_", df_region$gene[i], "_Coverage_plot.pdf"), width = 20, height = 20)
    plotTracks(list_to_plot, from = df_region$start[i], to = df_region$end[i], showTitle = TRUE)
    dev.off()
    
  }
  
}

##############################
# One plot = one condition
##############################

if (arguments$single == TRUE) {

  print("Chargement des fichiers")
  
  # Extract peak files directions            
  peak_file = readRDS(arguments$peak_file)
  name = str_extract(string = arguments$peak_file, pattern = "(?<=/genomic_ranges/).+(?=_threshold_10_ann.gr.rds)")
  bam_file = str_subset(list.files(arguments$bam_dir, pattern = name, full.names = TRUE), pattern = "-D[:digit:]{1}.bam(?!.bai|.nbreads|.qc)")
  
#########
# Tracking plots creation for each gene
#########

for(i in 1:nrow(df_region)) { # one plot per gene of the list

  print(paste("Boucle : gène", df_region$gene[i]))
  
  gr_region = GRanges(seqnames = df_region$chr[i],
                      ranges= IRanges(start = df_region$start[i],
                                      end = df_region$end[i]))

  # Chromosom with gene illustration
  # ********************************
  idxTrack <- IdeogramTrack(genome="hg19", chromosome= df_region$chr[i]) # drawing

  # Isoform tracking for the studied gene and promoters
  # ***************************************************
  ensGenes <- UcscTrack(genome = "hg19",
                        chromosome = df_region$chr[i],
                        background.title = "transparent",
                        alpha.title = 1,
                        col.title = "black",
                        cex.title = 1.5,
                        track = "ensGene",
                        from = df_region$start[i],
                        to = df_region$end[i],
                        trackType = "GeneRegionTrack",
                        rstarts = "exonStarts",
                        rends = "exonEnds",
                        transcript = "name",
                        strand = "strand",
                        fill = "lightgrey",
                        name = df_region$gene[i])

  # recupère seulement ce qui est commun entre les promoteurs et les régions définies dans grange
  prom_peak = subsetByOverlaps(prom_gene_fantom_gr, gr_region)

  df_prom_retained = tibble("start" = start(prom_peak),
                            "end" = end(prom_peak),
                            "chr" = as.character(seqnames(prom_peak)))

  if(nrow(df_prom_retained) != 0){  # faire apparaître les promoteurs sur le track des isoformes
    promoterTrack = HighlightTrack(trackList = ensGenes,
                                   start = df_prom_retained$start,
                                   end = df_prom_retained$end,
                                   chromosome = df_region$chr[i],
                                   col = "green",
                                   fill = "green")
  }

  print("Datatracking")
  
  # DataTracking
  dt = DataTrack(type = "histogram",
                 name = name,
                 background.title = "transparent",
                 fill.histogram = "#0072B2",
                 col.histogram = "#0072B2",
                 alpha.title = 1,
                 col.title = col_fun(name),
                 col.axis = "black",
                 cex.title = 1,
                 range = unlist(bam_file),
                 genome = "hg19",
                 ylim = c(0,df_region$max_read[i]),
                 window = -1,
                 chromosome = df_region$chr[i])


  # Create list of peaks detected in the gene region +/- 5kb
  peaks_in = subsetByOverlaps(peak_file, gr_region)

  print("plot")
  
  if (length(peaks_in) == 0){ # if not detected, no graph for this gene

    print("no_peaks_detected")

  }
  else # if detected, graph with highligted peaks
  {

    df_peaks_retained = tibble("start" = start(peaks_in),
                               "end" = end(peaks_in),
                               "chr" = as.character(seqnames(peaks_in)))

    ht = HighlightTrack(trackList = dt,
                        start = df_peaks_retained$start,
                        end = df_peaks_retained$end,
                        chromosome = df_region$chr[i])

    # Gather every tracking plot in a list
    if (exists("promoterTrack") == TRUE){
      list_to_plot = list(idxTrack, axTrack, blank_space, promoterTrack, blank_space, ht)
    } else{
      list_to_plot = list(idxTrack, axTrack, blank_space, ensGenes, blank_space, ht)
    }

    pdf(file =  paste0("D_Analysis/ATAC_gene_coverage/", name, "_", df_region$gene[i], "_Coverage_plot.pdf"), width = 8, height = 10)
    plotTracks(list_to_plot, from = df_region$start[i], to = df_region$end[i], showTitle = TRUE)
    dev.off()

  }

}
  
}

