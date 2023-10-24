#### Differential expression analysis ####

# When you see '## USER INPUT', this means you have to modify the code for your computer or dataset. All other code can be run as-is (i.e. you don't need to understand the code, just run it)

#### 2. Installing required packages ####

# **NOTE: this section only needs to be run once (or occasionally to update the packages)
# Install devtools
install.packages("devtools")
# Install R packages. This only needs to be run once.
bioconductor_packages <- c("DESeq2", "EnhancedVolcano", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Dr.eg.db", "org.Dm.eg.db")
cran_packages <- c("ggrepel", "ggplot2", "plyr", "reshape2", "readxl", "FactoMineR", "factoextra", "pheatmap")
# Compares installed packages to above packages and returns a vector of missing packages
new_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
# Install missing bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(new_packages)
# Install missing cran packages
if (length(new_cran_packages)) install.packages(new_cran_packages")
# Update all installed packages to the latest version
update.packages(bioconductor_packages, ask = FALSE)
update.packages(cran_packages, ask = FALSE")
