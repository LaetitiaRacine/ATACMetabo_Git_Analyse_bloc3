
# Test du package ChIPseeker avec un fichier Grange (MP_03h)
# http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#chip-peaks-coverage-plot

#**************************************************#
#*Chargement des librairies et fichiers de travail*#
#**************************************************#

# Chargement des librairies
# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# library(org.Hs.eg.db)
# library(stringr)


# Chargement des directions de travail
dir = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/"

# Chargement des peaks (format Grange) non annotés
# name_files = str_subset(list.files(path = paste0(dir, "A_Initial_data/merged_analysis/genomic_ranges"), 
#                                    pattern = ".gr.rds",
#                                    full.names = TRUE),
#                         pattern = "-D[:digit:].gr.rds")
# 
# files = list()
# for (i in 1:length(name_files)) {
#   name_list = str_extract(string = name_files[i], pattern = "(?<=ranges/).+(?=_D[:digit:])")
#   files = c(files, readRDS(name_files[i]))
#   names(files)[i] = name_list
# }
# 
# 
# peaks_gr = readRDS(paste0(dir, "A_Initial_data/merged_analysis/genomic_ranges/MP_03h_D1-D2-D3.gr.rds"))
# # peaks_gr_2 = readRDS(paste0(dir, "A_Initial_data/merged_analysis/genomic_ranges/MP_12h_D1-D2-D3.gr.rds"))
# 
# # Chargement des peaks (format Grange) pré-annotés
# peaks_gr_ann = readRDS(paste0(dir, "A_Initial_data/merged_analysis/genomic_ranges/MP_03h_D1-D2-D3_ann.gr.rds"))
# peaks_gr_ann


#**************************#
#* Fichiers en individuel *#
#**************************#

# Voir la distribution des annotations
# peaksAnno = annotatePeak(peaks_gr, tssRegion=c(-3000, 3000),
#                          TxDb=txdb, annoDb="org.Hs.eg.db")
# 
# plotAnnoPie(peaksAnno)
# plotAnnoBar(peaksAnno)
# vennpie(peaksAnno)
# 
# library(UpSetR)
# library(ggupset)
# library(ggplotify)
# # library(ggimage)
# upsetplot(peaksAnno)
# upsetplot(peaksAnno, vennpie = TRUE)


#**************************#
#* Comparaison de fichiers *#
#**************************#
# 
# promoter = getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
# tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
# plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
# plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
# tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
# 
# ## peak annotation comparision
# peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
#                        tssRegion=c(-3000, 3000), verbose=FALSE)
# 
# peakAnnoList_2DG <- lapply(files[1:3], annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
# peakAnnoList_2DG <- lapply(files[1:3], annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
# 
# plotAnnoBar(peakAnnoList_2DG)
# plotDistToTSS(peakAnnoList)
# 
# ## overlap of peaks and annotated genes
# genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# genes= lapply(peakAnnoList_2DG, function(i) as.data.frame(i)$geneId)
# 
# vennplot(genes)


#******************************************************************************#
#*** Upset plot without annotatePeak (ChIPseeker) *****************************#
#** https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html **#
#******************************************************************************#

library(dplyr)

# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# library(ComplexHeatmap)
# peaks_gr_ann = readRDS(paste0(dir, "A_Initial_data/merged_analysis/genomic_ranges/MP_03h_D1-D2-D3_ann.gr.rds"))
# peaks_ann = as.data.frame(peaks_gr_ann) %>% dplyr::select(-seqnames, -start, -end, -width, -strand)
# peaks_mat = t(as.matrix(peaks_ann))
# m = make_comb_mat(peaks_mat)
# pas possible trop de peaks


library(ggplot2)
library(tidyverse)

library(ggupset)   #https://github.com/const-ae/ggupset

peaks_gr_ann = readRDS(paste0(dir, "A_Initial_data/merged_analysis/genomic_ranges/MP_03h_D1-D2-D3_ann.gr.rds"))
peaks_ann = as.data.frame(peaks_gr_ann) %>% 
  dplyr::select(-seqnames, -start, -end, -width, -strand) %>%
  tibble::rownames_to_column(var = "peakID") 
peaks_pivot = pivot_longer(data = peaks_ann, cols = colnames(peaks_ann)[-1], names_to = "annotation", values_to = "value") %>%
  dplyr::filter(value == TRUE)
peaks_group_anno = peaks_pivot %>%
  group_by(peakID) %>%
  summarize(annotations = list(annotation)) 

plot1 = ggplot(peaks_group_anno, aes(x = annotations)) +
  geom_bar() +
  ggtitle("test") +
  ylab("Number of peaks") +
  scale_x_upset()
plot1

plot2 = ggplot(peaks_group_anno, aes(x = annotations)) +
  geom_bar() +
  ggtitle("test") +
  ylab("Number of peaks") +
  scale_x_upset(n_intersection = 50)
plot2



# UpSetR
a = peaks_pivot %>% 
  dplyr::select(-value) %>%
  unnest(cols = annotation) %>%
  mutate(annotMember=1) %>%
  pivot_wider(names_from = annotation, values_from = annotMember, values_fill = list(annotMember = 0)) %>%
  as.data.frame() %>%
  UpSetR::upset(sets = colnames(peaks_ann)[-1], 
                empty.intersections = NULL,
                nintersects = NA,
                order.by = "freq")
a

b = peaks_pivot %>% 
  dplyr::select(-value) %>%
  unnest(cols = annotation) %>%
  mutate(annotMember=1) %>%
  pivot_wider(names_from = annotation, values_from = annotMember, values_fill = list(annotMember = 0)) %>%
  as.data.frame() %>%
  UpSetR::upset(sets = colnames(peaks_ann)[-1], 
                empty.intersections = NULL,
                nintersects = 50,
                order.by = "freq") 
b
