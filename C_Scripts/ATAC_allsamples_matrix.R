
#**********************
# Libraries loading
#**********************

suppressPackageStartupMessages({
  suppressWarnings(library(Rsubread))
  suppressWarnings(library(dplyr))
})

dir = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/"
dir_initial_bam = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/A_Initial_data/merged_analysis/bam/"
dir_initial_gr = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/A_Initial_data/merged_analysis/genomic_ranges/"

#**********************
# Granges union for all samples 
#**********************

gr_list = lapply(paste0(dir_initial_gr, list.files(path = dir_initial_gr, pattern = "threshold_10_ann.gr.rds")), readRDS)
gr_intersection = Reduce(GenomicRanges::union, gr_list)
# saveRDS(gr_intersection, file=paste0(dir, "D_Analysis/ATAC_allsamples_matrix/union_gr_allsamples.gr.rds"))

df = data.frame(gr_intersection) %>% 
  tibble::rownames_to_column(var = "number") %>%
  dplyr::mutate(GeneID = paste0("peak_", number), .keep = "unused") %>%
  dplyr::rename(Chr = seqnames, Start = start, End = end, Strand = strand) %>%
  dplyr::select(GeneID, Chr, Start, End, Strand) 

#**********************
# readcount matrix for all sample union
#**********************

bam_files = str_subset(list.files(path = dir_initial_bam), pattern = "-D[:digit:].bam")
bam_files = paste0(dir_initial_bam, str_subset(bam_files, pattern = "bai|qc|nbreads", negate = TRUE))
  
# Count the number of reads, for each sample, present in the region resulting from the union
# readCount <- featureCounts(files = bam_files,
                           # annot.ext = df,  
                           # isPairedEnd = TRUE,
                           # nthreads = 1,
                           # countChimericFragments = FALSE,
                           # countMultiMappingReads = TRUE)

# saveRDS(readCount, file = paste0(dir, "D_Analysis/ATAC_allsamples_matrix/readcount_union_allsamples.rds"))

matrix_count = readCount$counts
# write.table(matrix_count, file = paste0(dir, "D_Analysis/ATAC_allsamples_matrix/featurecount_union_allsamples.txt"), sep = "\t", quote = FALSE)

count_df_read = tibble::rownames_to_column(data.frame(readCount$counts), "GeneID")
count_df_annotation = data.frame(readCount$annotation) %>% dplyr::select(-Strand, -Length)
count_df = left_join(count_df_annotation, count_df_read, by = "GeneID")
# write.table(count_df, file = paste0(dir, "D_Analysis/ATAC_allsamples_matrix/readcount_df_union_allsamples.csv"), sep = ";", row.names = FALSE) 

#**********************
# Analysis from : https://tobiasrausch.com/courses/atac/atac-seq-data-analysis.html#data-visualization
#**********************

# count_df = read.csv2(file = paste0(dir, "D_Analysis/ATAC_allsamples_matrix/readcount_df_union_allsamples.csv"))

library(DESeq2)
library(stringr)
library(tidyverse)

counts = count_df[,5:ncol(count_df)]

coldata = data.frame(file = names(counts), name = unlist(str_extract_all(names(counts), pattern = ".+(?=.bam)"))) %>%
  dplyr::mutate(donors = str_extract(name, pattern = "D[:digit:].+")) %>%
  dplyr::mutate(condtime = paste0(str_extract(name, pattern = ".+(?=h)"), "h")) %>% 
  tidyr::separate(col = condtime, into = c("condition", "time"), remove = TRUE) 
for (i in 1:nrow(coldata)) { if (coldata$condition[i] == "X2DG") coldata$condition[i] = "2DG" }
coldata = coldata %>%
  tidyr::unite(condition_time, condition, time, sep = "_", remove = FALSE) %>%
  remove_rownames %>%
  column_to_rownames(var = "file") %>%
  dplyr::select(-name)
  
dds = DESeqDataSetFromMatrix(
  countData = counts, 
  colData = coldata, 
  design = ~ condition)

dds = DESeq(dds)
cm = data.frame(counts(dds, normalized=TRUE))
rownames(cm) = paste0(count_df$Chr, '_', count_df$Start, '_', count_df$End)

pca = prcomp(t(cm))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
coldata = coldata %>% rownames_to_column(var = "sample")
pcaData=merge(pcaData, coldata)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

p=ggplot(data=pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = condition, shape = time), size=4) +
  scale_shape_manual(values=c(15, 16, 17, 7))
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
print(p)

q=ggplot(data=pcaData, aes(x = PC3, y = PC4)) + 
  geom_point(aes(colour = condition, shape = time), size=4) +
  scale_shape_manual(values=c(15, 16, 17, 7))
q=q+xlab(paste0("PC3: ", percentVar[3], "% variance"))
q=q+ylab(paste0("PC4: ", percentVar[4], "% variance"))
print(q)

varexp = data.frame(x=1:length(percentVar), y=percentVar)
varexp$x = factor(varexp$x)
ggplot(data=varexp, aes(x=x, y=y)) + geom_bar(stat="identity") + xlab("Principal Component") + ylab("Proportion of variation (%)")
