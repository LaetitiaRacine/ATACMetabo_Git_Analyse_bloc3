
# Gene coverage with Gviz
# https://www.bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html 


## Appel via docopt 

# arguments à appeler 
# arguments$prom = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/prom_gene_fantom_gr_corrected.rds"
# arguments$control = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/Xvivo_00h_D1-D2-D3_ann.gr.rds"
# arguments$bam_dir = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/bam/"
# arguments$gr_dir = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/"
# arguments$list_gene = c("GATA2","RUNX1","SMAD6","ERG","SPI1","CBFA2T3","FLI1","ZFPM1","HHEX","TAL1","GATA1")

#########
# Loading libraries
#########

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(GenomicRanges))
  suppressWarnings(library(Gviz))
})


#########
# Loading and storing files
#########

# Extract and organize gene information from genesymbol via Gviz
list_gene = c("GATA2","RUNX1","SMAD6","ERG","SPI1","CBFA2T3","FLI1","ZFPM1","HHEX","TAL1","GATA1")
max_read = c(25,20,15,10,30,15,20,20,40,40,10)

data(genesymbol, package = "biovizBase")
df_region = tibble()

for(i in 1:length(list_gene)) {
  region = genesymbol[list_gene[i]]
  region = range(region, ignore.strand = TRUE)
  region = region + 5000
  region = keepStandardChromosomes(region)
  df_region[i,"gene"] = list_gene[i]
  df_region[i,"chr"] = as.character(seqnames(region))
  df_region[i,"start"] = start(region)
  df_region[i,"end"] = end(region)
  df_region[i,"max_read"] = max_read[i]
}


prom_gene_fantom_gr = readRDS("/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/prom_gene_fantom_gr_corrected.rds")

# Extract peak control 00h directory and name and bam corresponding 
peak_control = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/Xvivo_00h_D1-D2-D3_ann.gr.rds"
peak_control_name = "Xvivo_00h_D1-D2-D3"
bam_control = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/bam/Xvivo_00h_D1-D2-D3.bam"

# Extract peak files directions            
peak_files = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/genomic_ranges/MP_03h_D1-D2-D3_ann.gr.rds"
peak_files_name = "MP_03h_D1-D2-D3"
bam_files = "/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/merged_analysis/bam/MP_03h_D1-D2-D3.bam"
# Extraire tous les Grange dans une liste si on veut tout représenter sur le même graph

# Load peak files in a list and name it
list_peaks_kept = lapply(X = c(peak_control, peak_files), FUN = readRDS)
names(list_peaks_kept) = c(peak_control_name, peak_files_name)
list_bam_files = list(bam_control, bam_files)
names(list_bam_files) = c(peak_control_name, peak_files_name)

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
      
      
      # DataTracking
        dt = DataTrack(type = "histogram",
                       name = paste0(names(list_peaks_kept[2]), "_track"),
                       background.title = "transparent",
                       fill.histogram = "#0072B2",
                       col.histogram = "#0072B2",
                       alpha.title = 1,
                       col.title = "black",
                       col.axis = "black",
                       cex.title = 1,
                       range = unlist(list_bam_files[2]),
                       genome = "hg19",
                       ylim = c(0,df_region$max_read[i]),
                       window = -1,
                       chromosome = df_region$chr[i])
      
      #########
      # 
      #########
        
      # Create list of peaks detected in the gene region +/- 5kb
        peaks_in = subsetByOverlaps(list_peaks_kept[[2]], gr_region)
        
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
          
          # Open a pdf file
          # pdf(file =  paste0("D_Analysis/ATAC_gene_coverage/", df_region$gene[i], "_Coverage_plot.pdf"), width = 8, height = 10)
          # Create a plot
          plotTracks(list_to_plot, from = df_region$start[i], to = df_region$end[i], showTitle = TRUE)
          # Close the pdf file
          # dev.off()
          
        }
        
    }

## Fonctionne pour une représentation mais il faut une boucle si on veut faire tout d'un coup !
 