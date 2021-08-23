
#########################
### Documentation
#########################

# Use of KaryoploteR from https://rdrr.io/bioc/karyoploteR/man/kpPlotRegions.html
# Other existing packages : ChIPseeker (covplot function) or ggbio

#########################
### Initialization : libraries loading, directories and function definitions
#########################

library(karyoploteR)
library(stringr)
library(tidyverse)

dir = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/"
dir_input = "/home/lracine/Documents/Git_Analyse_ATACMetabo_bloc3/A_Initial_data/merged_analysis/"

gr_list = str_subset(list.files(path = paste0(dir_input, "genomic_ranges"), pattern = "-D", full.names = T), pattern = "threshold_10_ann.gr.rds")

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
# One plot per sample with all genome
#########################

for (i in 1 : length(gr_list)) {
  
  gr = readRDS(gr_list[i])
  name_sample = str_extract(string = gr_list[i], pattern = "(?<=genomic_ranges/).+(?=_threshold)")
  
  # Peak distribution and peak density
  png(paste0(dir, "D_Analysis/ATAC_peaks_genome_distribution/plot_density_", name_sample, ".png"), width = 1000, height = 1200)
  kp = plotKaryotype(genome = "hg19", plot.type = 2) %>%
    kpPlotRegions(data=gr, avoid.overlapping = FALSE, col = col_fun(name_sample)) %>%
    kpAddMainTitle(main = paste0("Peaks distribution on Hg19 genome from ", name_sample)) %>%
    kpPlotDensity(data=gr, data.panel = 2)
  dev.off()
  
}

#########################
# One plot per chromosome with all samples
#########################

chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                "chr21", "chr22", "chrX", "chrY")

# Display settings : https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotParams/PlotParams.html 
pp = getDefaultPlotParams(plot.type=1)
pp$data1height = 10
pp$ideogramheight = 6
pp$bottommargin = 6
pp$data1inmargin = 6
  
for (i in 1:length(chromosomes)) {
  
  png(paste0(dir, "D_Analysis/ATAC_peaks_genome_distribution/plot_all_samples_distribution_density_",chromosomes[i], ".png"), width = 2000, height = 1000)
  kp = plotKaryotype(plot.type = 1, chromosomes=chromosomes[i], plot.params = pp)
  for (j in 1:length(gr_list)) {
      name_sample = str_extract(string = gr_list[j], pattern = "(?<=genomic_ranges/).+(?=_threshold)")
      kpPlotRegions(kp, data=readRDS(gr_list[j]), border=col_fun(name_sample), r0=j-1, r1=j-0.6)
      kpPlotDensity(kp, data=readRDS(gr_list[j]), data.panel = 1, r0=j-0.6, r1=j-0.2)
      kpAddLabels(kp, labels=name_sample, data.panel = 1, r0=j-1, r1=j-0.6) 
    }
  dev.off()
  
}

