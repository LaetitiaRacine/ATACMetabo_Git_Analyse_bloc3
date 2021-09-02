
#########################
### Command line to call the script in linux consol
#########################

"Create plot for gene coverage

Usage:
  ATAC_gene_biovizBase.R [options] <output_name>
  ATAC_gene_biovizBase.R -h | --help

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
  suppressWarnings(library(Gviz))
})


####################
## Script
####################

## Input
list_gene = c("GATA2","RUNX1","SMAD6","ERG","SPI1","CBFA2T3","FLI1","ZFPM1","HHEX","TAL1","GATA1")
max_read = c(25,20,15,10,30,15,20,20,40,40,10)

# Extract and organize gene information from genesymbol via Gviz
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

# Output
write.table(x = df_region, file = arguments$output_name, sep = ";")
            
            
            