
#########################
### Command line to call the script in linux consol
#########################

"Create plot by annotations types

Usage:
  ATAC_plot_peaks_annotation.R [options] <directory> <name_plot1> <name_plot2> <name_plot3> <name_plot4> <output.df_annot_wide> <output.df_annot_long>
  ATAC_plot_peaks_annotation.R -h | --help

Options:
  -h, --help              Show this screen

" -> doc

library(docopt)
arguments <- docopt(doc)

#########################
### Libraries loading
#########################

suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(stringr))
  #suppressWarnings(library(tidyverse))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(viridis))
})

############################################################
### Input loading and organizing in data.frame
############################################################

peaks = str_subset(list.files(path = arguments$directory), pattern = "threshold_[:digit:]{2,3}_ann.gr.rds")

df_peaks_annotation = data.frame()
df_peaks_recap = data.frame()

for (i in 1:length(peaks)) {
  temp = as.data.frame(readRDS(paste0(arguments$directory,peaks[i]))) %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>%
    dplyr::mutate(sample = str_extract(string = peaks[i], pattern = "[:graph:]+(?=_threshold)"))
  temp2 = temp %>% 
    group_by(sample) %>%
    summarise(across(.cols = 1:(ncol(temp)-1), .fns = sum)) %>%
    tidyr::separate(col = sample, into = c("condition","time","donor"), sep = "_", remove = TRUE) %>%
    dplyr::select(-donor)
  df_peaks_annotation = rbind(df_peaks_annotation, temp2)
  temp3 = temp2 %>%
    dplyr::select(condition, time) %>%
    dplyr::mutate(nbpeaks = nrow(temp))
  df_peaks_recap = rbind(df_peaks_recap, temp3)
}

df_peaks_annotation_l = pivot_longer(data = df_peaks_annotation, cols = 3:ncol(df_peaks_annotation), names_to = "annotation", values_to = "count")

#####################
### Plots creation 
#####################

nbpeaks_pertime = ggplot(df_peaks_recap,aes(x = time,  y = nbpeaks, label = nbpeaks, fill = condition)) +
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
nbpeaks_pertime

nbpeaks_percond = ggplot(df_peaks_recap,aes(x = condition,  y = nbpeaks, label = nbpeaks, fill = time)) +
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
  scale_fill_viridis_d()
nbpeaks_percond

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
  scale_fill_viridis_d()
annot_percond

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
annot_pertime

#####################
### Outputs' saving 
#####################

ggsave(plot = nbpeaks_pertime,
       filename = arguments$name_plot1,
       width = 16*0.75, height = 9*0.75)

ggsave(plot = nbpeaks_percond, 
       filename = arguments$name_plot2,
       width = 16*0.75, height = 9*0.75)

ggsave(plot = annot_pertime, 
       filename = arguments$name_plot3,
       width = 16*0.75, height = 9*0.75)

ggsave(plot = annot_percond, 
       filename = arguments$name_plot4,
       width = 16*0.75, height = 9*0.75)

write.table(df_peaks_annotation, file = arguments$output.df_annot_wide, sep = ";", row.names = FALSE) 
write.table(df_peaks_annotation_l, file = arguments$output.df_annot_long, sep = ";", row.names = FALSE) 
