#**********************
# Command line to call the script in linux consol
#**********************

"Create plot by annotations types

Usage:
  Plot_Annotation.R [options] <directory>
  Plot_Annotation.R -h | --help

Options:
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
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(stringr))
})


#*********************
#
#*********************

# arguments$directory
test = "/home/lracine/Documents/Bloc2_data_ATACMetabo/D_Analysis_merged/genomic_ranges/"
# travail avec les fichiers threshold 10 pour tester 

# Load files of interest
peak_files = dir(path = test, pattern = "10_ann.gr.rds", full.names = TRUE)
peak_list = lapply(X = peak_files, FUN = readRDS)
condition_names = dir(path = test, pattern = "10_ann.gr.rds")
condition_names = str_extract(condition_names, pattern = ".+(?=_threshold)")
names(peak_list) = condition_names 


peak_list_static = peak_list[str_detect(names(peak_list), pattern = "intersection")]





##########
# 2 : Count peaks occurences in annotation types
#########

# Define annotation types
annotations_types = names(mcols(peak_list[[1]]))# 1 or whatever, they're all the same
count_all_df = tibble()

for(i in 1:length(peak_list)){
  
  #extract sum of TRUE by annotation type
  count_ann = as.data.frame(mcols(peak_list[[i]])) %>%
    summarise_all(list(sum = sum)) %>% #Sum all true values in each colon (result = data frame 1 row and n columns)
    t() %>% # transpose (devient matrix)
    as.data.frame() %>%
    mutate(condition = names(peak_list)[i] , annotation = annotations_types)
  
  
  ### create df combining sum of TRUE for each annotation type by time
  count_all_df = rbind(count_all_df, count_ann)
  
}

# Clean the data frame and reorder it
count_all_df = count_all_df[,c(3,2,1)]
colnames(count_all_df)[3] = "Value"

#!allow to separate static and differential analysis
count_all_df = count_all_df %>%
  mutate(analysis = ifelse(str_detect(condition,pattern = "vs"), "differential", "static"))

#Create a dataframe
count_all_df_static = count_all_df %>%
  filter(analysis == "static") %>%
  mutate(new_condition = str_sub(condition, 20))

count_all_df_differential = count_all_df %>%
  filter(analysis == "differential") %>%
  mutate(new_condition = str_sub(condition, 7))

list_count_all = list("count_all_df_static" = count_all_df_static, "count_all_df_differential" =  count_all_df_differential)

