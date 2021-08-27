
#################
# Documentation #
#################

# The original "prom_gene_fantom_gr_rdata" send from Ravi (I think) lack HHEX FANTOM5 information.
# We create a corrected version with the missing datas.

###########################
# Libraries and functions #
###########################

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

###########################
# Script #
###########################

# Create a Grange with HHEX information
HHEX_prom = GRanges(seqnames = "chr10",
                    ranges = IRanges(start = c(94449649,94449675,94449703,94451574), end = c(94449664,94449694,94449718,94451587)), 
                    strand = "*")
mcols(HHEX_prom) = tibble(gene = "HHEX")

# Load prom_gene_fantom_gr and add HHEX gr 
prom_gene_fantom_gr = loadRData("/home/lracine/Documents/ATACMetabo_Git_Analyse_bloc3/A_Initial_data/prom_gene_fantom_gr.rdata")
prom_gene_fantom_gr_corrected = c(prom_gene_fantom_gr, HHEX_prom)
prom_gene_fantom_corrected = sort(prom_gene_fantom_gr_corrected)

# Save the new file
saveRDS(prom_gene_fantom_corrected, file = "A_Initial_data/prom_gene_fantom_gr_corrected.rds")
