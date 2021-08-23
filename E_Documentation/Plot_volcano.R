#####################
##### Libraries #####
#####################

library(GenomicRanges)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(ggrepel)


######################
##### Functions #####
######################

#Loads an RData file, and assign it a new name
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

######################
####### Script #######
######################

#the_dir = ""/home/doctorant/Documents/Temp_ATAC_0620/"
the_dir = "~/Bureau/ATAC_2020/"
setwd(the_dir)

########### 
# 1 : Load and clean files of interest  
###########

list_d2res_files = dir(path = paste0(the_dir,"results/DEseq_output"))

for( i in 1:length(list_d2res_files)){

d2res = read.delim(paste0(the_dir,"results/DEseq_output/",list_d2res_files[[i]]))

condition = str_sub(list_d2res_files[[i]], 16, -5)
list_peak_files = dir(path = paste0(the_dir,"results/genomic_ranges/annotated_genomic_ranges"), pattern = condition, full.names = TRUE)

# Make a list with both  Grange correpsonding to the condition (inc and dec)
DEseq_peaks_list = lapply(list_peak_files, loadRData)

# Concatenate those 2 Grange
all_deseq_peaks = c(DEseq_peaks_list[[1]], DEseq_peaks_list[[2]])

###########
# 2 : Adding peaks features information in d2re
###########

intergenic_peaks = all_deseq_peaks[all_deseq_peaks$Intergenic == TRUE]
intergenic_peaks = rownames(mcols(intergenic_peaks))

TSS_peaks = all_deseq_peaks[all_deseq_peaks$TSS_mp1kb == TRUE]
TSS_peaks = rownames(mcols(TSS_peaks))

CTCF_peaks = all_deseq_peaks[all_deseq_peaks$CTCF == TRUE]
CTCF_peaks = rownames(mcols(CTCF_peaks))

#Add a column "region" in d2res where Intergenic, TSS and  CTCF are notified
d2res$region = ""
d2res[d2res$name %in% intergenic_peaks,]$region = "Intergenic"
d2res[d2res$name %in%CTCF_peaks,]$region = "CTCF"
d2res[d2res$name %in%TSS_peaks,]$region = "TSS_mp1kb"


###########
# 3 : Statistics about dynamic
###########

plotDF<- d2res[,c(1,6,9,11,12)]

regulation_statistics <- table(plotDF$regulation)

total = regulation_statistics[1] + regulation_statistics[2]
per_Open = round(regulation_statistics[2] /total,3)*100
per_Close = round(regulation_statistics[1]/total,3)*100

###########
# 4 : Focusing on gene of inereest
###########

# On charge la base de données de promoter que Ravi nous a fourni
prom_gene_fantom_gr = loadRData(paste0(the_dir,"data/prom_gene_fantom_gr.rdata"))

# Adding FANTOM5 information for HHEX (lacking in original dataset, don't know why ?)

HHEX_prom = GRanges(seqnames = "chr10",
                    ranges = IRanges(start = c(94449649,94449675,94449703,94451574), end = c(94449664,94449694,94449718,94451587)), 
                    strand = "*")

mcols(HHEX_prom) = tibble(gene = "HHEX")
prom_gene_fantom_gr = c(prom_gene_fantom_gr, HHEX_prom)
prom_gene_fantom_gr = sort(prom_gene_fantom_gr)

# List of genes of interest (want to see them appear on the plot)
list_gene = c("GATA2","GATA1","RUNX1","SMAD6","ERG","SPI1","CBFA2T3","FLI1","ZFPM1","HHEX","TAL1")

# Retrieve all the promoters for the gene of interest
prom_gene_fantom_gr = prom_gene_fantom_gr[prom_gene_fantom_gr$gene %in% list_gene]

# Tranform the D2res table in a Grange for overlapping calculation
volcano_genes_gr = GRanges(seqnames = d2res$Chr,
  ranges = IRanges(start = d2res$Start, end = d2res$End))

# Retrieves regions overlapping between D2res regions and promoters for gene of interest
peaks_hemato_TSS = findOverlaps(volcano_genes_gr, prom_gene_fantom_gr)
overlap_peaks_nb = queryHits(peaks_hemato_TSS)
overlap_prom_nb = subjectHits(peaks_hemato_TSS)

cor_peak_prom_gene = tibble("peaks_nb" = plotDF[overlap_peaks_nb,]$name,
                            "gene_name" = prom_gene_fantom_gr$gene[overlap_prom_nb])

cor_peak_prom_gene = distinct(cor_peak_prom_gene)


### Traçage du volcanoplot change-nochange
downColor <- "red"
upColor <- "green4"
selected_color <- "black"

# Redefine levels
plotDF$regulation <- factor(plotDF$regulation, levels = c("gained-close", "gained-open", "selected")) # All levels possible (selected used later)
plotDF$region = factor(plotDF$region, levels = c("TSS_mp1kb","Intergenic","","CTCF"))

plotDF$Gene = "" # Empty column creation for gene name associated to the peak
overlap_peaks = which(plotDF$name %in% cor_peak_prom_gene$peaks_nb)
plotDF[overlap_peaks,]$Gene = cor_peak_prom_gene$gene_name
plotDF[overlap_peaks,]$regulation = "selected"


plotSubTitle <- paste(paste0("gained-closed regions = ",regulation_statistics[1]," (",per_Close,"%)"),
                      paste0("gained-open regions = ",regulation_statistics[2]," (",per_Open,"%)"),
                      paste0(condition," (total = ", total, ")"),
                      sep = "      ")
xlimit = 7

p <- ggplot(plotDF, aes(x = log2FoldChange, y = -1 * log10(pvalue), col=regulation, label = Gene)) +
  geom_point(size = 2, alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
  scale_color_manual(name = "Regulation", values = c(downColor,upColor,selected_color)) +
  scale_shape_manual(name = "Regulation" , values = c(20,20,20)) +
  scale_x_continuous(limits = c(-xlimit, xlimit)) +
  ylim(NA,20) +
  labs(subtitle = plotSubTitle,x = "log2(FoldChange)", y = "-log10(Pvalue)") +
  theme_tufte()+
  theme(axis.line.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    axis.title = element_text(size = 10),
    plot.subtitle = element_text(size = 11,hjust = 0.5)) +
  geom_text_repel(data = subset(plotDF, Gene != ""),
    # color = "transparent",
    box.padding = 3,
    seed = 9,
    fontface = 'bold',
    segment.color = '#336699',
    force = 9)

  # print(p)
      ggsave(p,
            units = "mm",
            width = 170, 
            dpi = 600,
            filename = paste0(the_dir, "results/plot/volcano_", condition, ".jpeg"))


  ### Only plot TSS and Intergenic regions

plotDF_2 = plotDF %>%
          dplyr::filter(region =="TSS_mp1kb" | region == "Intergenic")

p <- ggplot() +
    geom_point(data = plotDF_2,
      aes(x = log2FoldChange, y = -1 * log10(pvalue),
      color = region),
      alpha = 0.6,
      size = 3)+
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
    geom_text_repel(data = subset(plotDF_2, Gene %in% list_gene),
        aes(x = log2FoldChange, y = -1 * log10(pvalue),
          color = region,
          label = Gene),
        box.padding = 3,
        color = "transparent",
        seed = 9,
        fontface = 'bold',
        segment.color = '#909498',
        force = 9)+
    scale_color_manual(name = "Genomic feature", values = c("#0072B2","#F0E442","white","white")) +
    scale_x_continuous(limits = c(-xlimit, xlimit)) +
    ylim(NA,20) +
  labs(subtitle = paste0("pValue < 0.05 : ", nrow(plotDF[plotDF$pvalue < 0.05,]), " peaks"),
    x = "log2(FoldChange)",
    y = "-log10(Pvalue)") +    theme_tufte()+
    theme(axis.line.y = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      axis.title = element_text(size = 10),
      plot.subtitle = element_text(size = 11,hjust = 0.5))
  
 # print(p)
      ggsave(p,
            units = "mm",
            width = 170, 
            dpi = 600,
            filename = paste0(the_dir, "results/plot/volcano_TSS_Intergenic_", condition, ".jpeg"))
      
      p <- ggplot() +
    geom_point(data = plotDF_2,
      aes(x = log2FoldChange, y = -1 * log10(pvalue),
      color = region),
      alpha = 0.6,
      size = 3)+
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed")+
    scale_color_manual(name = "Genomic feature", values = c("#0072B2","#F0E442","white","white")) +
    scale_x_continuous(limits = c(-xlimit, xlimit)) +
    ylim(NA,20) +
  labs(subtitle = paste0("pValue < 0.05 : ", nrow(plotDF[plotDF$pvalue < 0.05,]), " peaks"),
    x = "log2(FoldChange)",
    y = "-log10(Pvalue)") +    theme_tufte()+
    theme(axis.line.y = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      axis.title = element_text(size = 10),
      plot.subtitle = element_text(size = 11,hjust = 0.5))
  
 # print(p)
      ggsave(p,
            units = "mm",
            width = 170, 
            dpi = 600,
            filename = paste0(the_dir, "results/plot/volcano_TSS_Intergenic_", condition, "without_lines.jpeg"))
}
