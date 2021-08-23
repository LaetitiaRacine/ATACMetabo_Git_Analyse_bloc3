
## Idée : ajouter des options d'appel dans la rule (voir report_merged_plot bloc2) et ajouter des if dans le code 


#########################
### Command line to call the script in linux consol
#########################

"Create plot by annotations types

Usage:
  ATAC_plot_peaks_annotation.R <type> <gr_info> <output> 
  ATAC_plot_peaks_annotation.R <type> <df_peaks> <df_annot> <facet> <plot_name> 
  ATAC_plot_peaks_annotation.R -h | --help

Options:
  -h, --help              Show this screen

" -> doc

library(docopt)
arguments <- docopt(doc)

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

#########################
### Libraries loading
#########################

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(stringr))
  suppressWarnings(library(tidyverse))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(viridis))
  suppressWarnings(library(ggupset))
  suppressWarnings(library(grid))
})

############################################################
### Table : Nb peaks per annotation for each sample
############################################################

if (arguments$type == "table") {
  
  peaks = str_subset(list.files(path = arguments$gr_info), pattern = "threshold_[:digit:]{2,3}_ann.gr.rds")
  
  df_peaks_annotation = data.frame()

  for (i in 1:length(peaks)) {
    temp = as.data.frame(readRDS(paste0(arguments$gr_info, peaks[i]))) %>%
      dplyr::select(-seqnames, -start, -end, -width, -strand) %>%
      dplyr::mutate(sample = str_extract(string = peaks[i], pattern = "[:graph:]+(?=_threshold)"))
    temp2 = temp %>% 
      group_by(sample) %>%
      summarise(across(.cols = 1:(ncol(temp)-1), .fns = sum)) %>%
      tidyr::separate(col = sample, into = c("condition","time","donor"), sep = "_", remove = TRUE)
    df_peaks_annotation = rbind(df_peaks_annotation, temp2)
  }
  
  write.table(df_peaks_annotation, file = arguments$output, sep = ";", row.names = FALSE) 

}

#####################
### Plot : Nb peaks per time or condition
#####################

if (arguments$type == "nbpeaks") {
  
  df_peaks = read.csv2(arguments$df_peaks, sep = ";") %>% dplyr::filter(threshold == "10")
  
  if (arguments$facet == "pertime") {
    
    nbpeaks_pertime = ggplot(df_peaks,aes(x = time,  y = nbpeaks, label = nbpeaks, fill = condition)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
      geom_text(position = position_dodge2(width = 0.9, preserve = "single"), angle = 90, hjust=1.5) +
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
            strip.text.x = element_text(size=11, color="black", face="bold.italic"),
            axis.line.y = element_line(colour = "black"),
            axis.text.y = element_text(size = 11, colour = "black"),
            axis.ticks.y = element_line() ,
            axis.title.y = element_text(vjust = 1 ,size = 14),
            axis.title.x = element_blank(),
            axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black"),
            axis.ticks = element_blank()) +
      labs(y = "Number of peaks detected")
    
    ggsave(plot = nbpeaks_pertime,
           filename = arguments$plot_name,
           width = 16*0.75, height = 9*0.75)
    
  }

  if (arguments$facet == "percond") {
    
    nbpeaks_percond = ggplot(df_peaks,aes(x = condition,  y = nbpeaks, label = nbpeaks, fill = time)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
      geom_text(position = position_dodge2(width = 0.9, preserve = "single"), angle = 90, hjust=1.5) +
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
            strip.text.x = element_text(size=11, color="black", face="bold.italic"),
            axis.line.y = element_line(colour = "black"),
            axis.text.y = element_text(size = 11, colour = "black"),
            axis.ticks.y = element_line() ,
            axis.title.y = element_text(vjust = 1 ,size = 14),
            axis.title.x = element_blank(),
            axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black"),
            axis.ticks = element_blank()) +
      labs(y = "Number of peaks detected") +
      scale_fill_viridis_d(alpha = 0.7)

    ggsave(plot = nbpeaks_percond, 
           filename = arguments$plot_name,
           width = 16*0.75, height = 9*0.75)
    
  }

}

#####################
### Plot : Nb peaks per time or condition - facet by annotations
#####################

if (arguments$type == "nbpeaks_annotation") {
  
  df_peaks_annotation = read.csv2(arguments$df_annot, sep = ";")
  df_peaks_annotation_l = pivot_longer(data = df_peaks_annotation, cols = 4:ncol(df_peaks_annotation), names_to = "annotation", values_to = "count")
  
  if (arguments$facet == "pertime") {
    
    annot_pertime = ggplot(df_peaks_annotation_l,aes(x = time,  y = count, fill = condition)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
            strip.text.x = element_text(size=11, color="black", face="bold.italic"),
            axis.line.y = element_line(colour = "black"),
            axis.text.y = element_text(size = 11, colour = "black"),
            axis.ticks.y = element_line() ,
            axis.title.y = element_text(vjust = 1 ,size = 14),
            axis.title.x = element_blank(),
            axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black")) +
      facet_wrap(annotation~., scale = "free_x")
    
    ggsave(plot = annot_pertime, 
           filename = arguments$plot_name,
           width = 16*0.75, height = 9*0.75)
    
  }
  
  if (arguments$facet == "percond") {
    
    annot_percond = ggplot(df_peaks_annotation_l,aes(x = condition,  y = count, fill = time)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5, face = "bold", colour= "black", size = 16),
            strip.text.x = element_text(size=11, color="black", face="bold.italic"),
            axis.line.y = element_line(colour = "black"),
            axis.text.y = element_text(size = 11, colour = "black"),
            axis.ticks.y = element_line() ,
            axis.title.y = element_text(vjust = 1 ,size = 14),
            axis.title.x = element_blank(),
            axis.text.x = element_text(vjust = -0.5, size = 11, colour = "black")) +
      facet_wrap(annotation~., scale = "free_x") +
      scale_fill_viridis_d(alpha = 0.7)

    ggsave(plot = annot_percond, 
           filename = arguments$plot_name,
           width = 16*0.75, height = 9*0.75)
  }
  
}

#####################
### Plot : Upset plot, annotations and intersection
#####################

if (arguments$type == "upset" | arguments$type == "upset_zoom") {

  peaks_ann = as.data.frame(readRDS(arguments$gr_info)) %>% 
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>%
    tibble::rownames_to_column(var = "peakID") 
  peaks_pivot = pivot_longer(data = peaks_ann, cols = colnames(peaks_ann)[-1], names_to = "annotation", values_to = "value") %>%
    dplyr::filter(value == TRUE)
  peaks_group_anno = peaks_pivot %>%
    group_by(peakID) %>%
    summarize(annotations = list(annotation)) 
  name_sample = str_extract(arguments$gr_info, "(?<=ranges/).+(?=_ann)")
    
  if (arguments$type == "upset") {
    
    plot = ggplot(peaks_group_anno, aes(x = annotations)) +
      geom_bar(fill = col_fun(arguments$gr_info), color = col_fun(arguments$gr_info)) +
      ggtitle(paste0("Peaks distribution in genomic features - ", name_sample)) +
      ylab("Number of peaks") +
      scale_x_upset()
    ggsave(plot = plot, filename = arguments$output, width = 16*0.75, height = 9*0.75)
    
    plot_name = arguments$output
    str_sub(plot_name, str_locate(plot_name, "plot_")[2]+1, str_locate(plot_name, "plot_")[2]) <- "set_size_"
    png(plot_name, width = 2000, height = 1000)
    print(peaks_pivot %>%
            dplyr::select(-value) %>%
            unnest(cols = annotation) %>%
            mutate(annotMember=1) %>%
            pivot_wider(names_from = annotation, values_from = annotMember, values_fill = list(annotMember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(sets = colnames(peaks_ann)[-1],
                          empty.intersections = NULL,
                          nintersects = NA,
                          order.by = "freq",
                          main.bar.color = col_fun(arguments$gr_info))
    )
    grid.text(paste0("Peaks distribution in genomic features - ", name_sample), x = 0.65, y=0.95, gp=gpar(fontsize=20))
    dev.off()

  }
  
  
  if (arguments$type == "upset_zoom") {
    
    plot = ggplot(peaks_group_anno, aes(x = annotations)) +
      geom_bar(fill = col_fun(arguments$gr_info), color = col_fun(arguments$gr_info)) +
      ggtitle(paste0("Peaks distribution in genomic features - ", name_sample)) +
      ylab("Number of peaks") +
      scale_x_upset(n_intersection = 50)
    ggsave(plot = plot, filename = arguments$output, width = 16*0.75, height = 9*0.75)

    plot_name = arguments$output
    str_sub(plot_name, str_locate(plot_name, "plot_")[2]+1, str_locate(plot_name, "plot_")[2]) <- "set_size_"
    png(plot_name, width = 2000, height = 1000)
    print(peaks_pivot %>%
            dplyr::select(-value) %>%
            unnest(cols = annotation) %>%
            mutate(annotMember=1) %>%
            pivot_wider(names_from = annotation, values_from = annotMember, values_fill = list(annotMember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(sets = colnames(peaks_ann)[-1],
                          empty.intersections = NULL,
                          nintersects = 50,
                          order.by = "freq",
                          main.bar.color = col_fun(arguments$gr_info))
    )
    grid.text(paste0("Peaks distribution in genomic features - ", name_sample), x = 0.65, y=0.95, gp=gpar(fontsize=20))
    dev.off()
    
  }
  
  
}
