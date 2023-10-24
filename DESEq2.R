#### 3. Loading required packages ####

# This section needs to be run every time
# Load packages
bioconductor_packages <- c("DESeq2", "EnhancedVolcano", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Dr.eg.db", "org.Dm.eg.db")
cran_packages <- c("ggrepel", "ggplot2", "plyr", "reshape2", "readxl", "FactoMineR", "factoextra", "pheatmap")
lapply(cran_packages, require, character.only = TRUE)
lapply(bioconductor_packages, require, character.only = TRUE)

#### 4. Import your count data ####

# Make sure you have: a) your count table (salmon.merged.gene_counts.tsv file, if you used Nextflow nfcore/rnaseq to analyse your data). Copy this to a subdirectory called 'data'. b) your metadata file. This should be either an Excel file called 'metadata.xlsx' or a tab-separated text file called 'metadata.txt'. It needs 3 columns called 'sample_name', 'sample_ID' and 'group'. The sample names should be EXACTLY the same as the names in the count table. These names are often uninformative and long, so the 'sample_ID' is the sample labels you want to put on your plots. E.g. if you have a 'high fat' group, you might want to rename the samples HF1, HF2, HF3, etc)

## USER INPUT
# Set working directory. 
# Change this to your working directory (In the RStudio menu: Session -> Set working directory -> Choose working directory)
setwd("H:/workshop/RNAseq")

# Import your count data. make sure you've created a 'data' subdirectory and put the count table file there.
metacountdata <- read.table("./data/salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1)

# Import metadata. Again, need a metadata.xlsx file in the data subdirectory.
meta <- read_excel("./data/metadata.xlsx")

# Remove 1st columns of metadata (gene_name)
counts <- metacountdata[ ,2:ncol(metacountdata)]

# Rename sample names to new sample IDs
counts <- counts[as.character(meta$sample_name)]
colnames(counts) <- meta$sample_ID

# Counts need to be rounded to integers
counts <- ceiling(counts)

#### 5. Outliers and batch effects ####

# This section normalises and transforms the count data so that it can be plotted on a PCA plot and a heatmap

## USER INPUT
# Choose the groups you want to plot in a PCA/Heatmap. You can select any 2 or more of the groups (or all of the groups) you have in your 'groups' column of your metadata table
# To see what groups are present, run the following:
unique(meta$group)
# Now add which groups you want to plot (i.e. replace the groupnames below, and add more, separated by a comma and in "quotes", as needed). NOTE: R is case-sensitive, so these group names must be named EXACTLY the same as in the metadata table.
plotgroups <- c("Differentiated_cells", "Basal_cells")

# Pull out only the counts from the above groups
groupcounts <- counts[meta$group %in% plotgroups]

# Normalise counts by library size, using DeSeq2's estimateSizeFactors() function. Note that DeSeq2 does this internally during DEG calling. The normalisation below is done separately for PCA and density plotting.
# Set up the initial DeSeq2 experimental parameters.
condition <- factor(1:length(groupcounts))
# Set up the column data. A data frame of sample ID's and conditions
coldata <- data.frame(row.names=colnames(groupcounts), condition)
# Set up the DeSeq2 data set structure
f <- DESeqDataSetFromMatrix(countData = groupcounts, colData = coldata, design= ~ condition)
# Estimate the size factors. See DeSeq2 manual for details
f <- estimateSizeFactors(f)
# Size factors can be viewed by: sizeFactors(f)

# Multiply each row (sample) by the corresponding size factor
subcount_norm <- as.matrix(groupcounts) %*% diag(sizeFactors(f))
# Re-add column names
colnames(subcount_norm) <- colnames(groupcounts)

## Remove low coverage transcripts (mean count < 10) ##

# Find the mean of each row (and output as a data frame object)
means <- as.data.frame(rowMeans(subcount_norm))
# Then join the means data with the counts
means <- cbind(means, subcount_norm)
# Then subset out only genes with mean > 10
data <- subset(means, means[ , 1] > 10)
# Remove the means column
data <- data[,-1]

# Transform data
data_log <- vst(round(as.matrix(data)))
# Transformation can create some infinite values. Can't generate PCA data on these. Can see how many by: sum(sapply(data_log, is.infinite))
# To remove infinite rows, use 'is.finte' or '!is.infinite'
data_log <- data_log[is.finite(rowSums(data_log)),]
colnames(data_log) <- colnames(groupcounts)

### Set up the PCA plot base data ###

# We're using the FactoMineR package to generate PCA plots (http://factominer.free.fr/index.html)

# Need to transpose the data first
data_log_t <- t(data_log)
# Add the group data
data_log_t_vars <- data.frame(meta$group[meta$group %in% plotgroups], data_log_t)
# Generate the PCA data using FactoMineR package
res.pca <- PCA(data_log_t_vars, quali.sup = 1, graph=FALSE)

## Set up the dendogram/heatmaps base data ##

# Calculate the distance matrix:
distance_matrix <- as.matrix(dist(t(data_log)))


#### 5a. PCA plot ####

# Generate the PCA plot. Groups are shaded with ellipses at 95% confidence level. NOTE: at least 4 replicates need to be in a group for an ellipses to be drawn.
# NOTE: change the group point colours by changing 'palette = ' below. Use the 'RColourBrewer' colour names (https://r-graph-gallery.com/38-rcolorbrewers-palettes.html). For example, if you are plotting 3 groups and choose palette = "Set1", this will use the first 3 colours from the Set1 colour palette.
p <- fviz_pca_ind(res.pca,
                  geom.ind = "point", # show points only (but not "text")
                  col.ind = meta$group[meta$group %in% plotgroups], # color by groups
                  pointsize = 5, title = "", legend.title = "Treatment groups", palette = "Dark2",
                  addEllipses = TRUE, ellipse.type = "t", ellipse.level = 0.95) + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title=element_text(size=16), axis.text=element_text(size=14))

p

# Output as publication quality (300dpi) tiff and pdf.
# This will name your output files with the treatment groups you selected.

# Create a 'figures' subdirectory where all figures will be output
dir.create("figures", showWarnings = FALSE)
# Create a (300dpi) tiff
ggsave(file = paste0("./figures/PCA_", paste(plotgroups, collapse = "_Vs_"), ".tiff"), dpi = 300, compression = "lzw", device = "tiff", width = 10, height = 8, plot = p)
# Create a pdf
ggsave(file = paste0("./figures/PCA_", paste(plotgroups, collapse = "_Vs_"), ".pdf"), device = "pdf", width = 10, height = 8, plot = p)


#### 5b. Samples heatmap and dendrogram ####

# This section plots a heatmap and dendrogram of pairwise relationships between samples. In this way you can see if samples cluster by treatment group.

# See here: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

# Define annotation column
annot_columns <- data.frame(meta$group[meta$group %in% plotgroups])
# Make the row names the sample IDs
row.names(annot_columns) <- meta$sample_ID[meta$group %in% plotgroups]
colnames(annot_columns) <- "Treatment groups"
# Need to factorise it
annot_columns[[1]] <- factor(annot_columns[[1]])

# Generate dendrogram and heatmap
pheatmap(distance_matrix, color=colorRampPalette(c("white", "#9999FF", "#990000"))(50), cluster_rows = TRUE, show_rownames = TRUE, treeheight_row = 0, treeheight_col = 70, fontsize_col = 12, annotation_names_col = F, annotation_col = annot_columns, filename = paste0("./figures/Pairwise_sample_heatmap_", paste(plotgroups, collapse = "_Vs_"), ".tiff"))

# Notes about heatmap colours.
# You can change the colours used in the heatmap itself by changing the colour names (color=colorRampPalette....)
# If you want to change the annotation colours, see here: https://zhiganglu.com/post/pheatmap_change_annotation_colors/


#### 6. Differential expression analysis ####

# In this section we use the Deseq2 package to identify differentially expressed genes.

## USER INPUT
# Choose the treatment groups you want to compare. 
# To see what groups are present, run the following:
unique(meta$group)
# Enter which groups you want to compare (two groups only). BASELINE OR CONTROL GROUP SHOULD BE LISTED FIRST.
degroups <- c("Basal_cells", "Differentiated_cells")

# From the count table, pull out only the counts from the above groups
expdata <- as.matrix(counts[,meta$group %in% degroups])

# Set up the experimental condition
# 'factor' sets up the reference level, i.e. which is the baseline group (otherwise the default baseline level is in alphabetic order)
condition <- factor(meta$group[meta$group %in% degroups], levels = degroups)
# Type 'condition' in the console to see is the levels are set correctly

# Set up column data (treatment groups and sample ID)
coldata <- data.frame(row.names=colnames(expdata), condition)

# Create the DESeq2 dataset (dds)
dds <- DESeq2::DESeqDataSetFromMatrix(countData=expdata, colData=coldata, design=~condition)
dds$condition <- factor(dds$condition, levels = degroups)

# Run DESeq2 to identify differentially expressed genes
deseq <- DESeq(dds)

# Extract a results table from the DESeq analysis
res <- results(deseq)

# Reorder results by adjusted p vales, so that the most signififcantly DE genes are at the top
res <- res[order(res$padj), ]

# You can do a summary of the results to see how many significantly (alpha=0.05, adjust to 0.01 if needed) upregulated and downregulated DE genes were found
summary(res, alpha=0.05)

# Convert from DESeq object to a data frame.
res <- data.frame(res)

# Look at the top 6 DE genes
head(res)


#### 6a. Annotating your DE genes ####

# Annotation packages for human (org.Hs.eg.db), mouse (org.Mm.eg.db), rat (org.Rn.eg.db), E. coli strain K12 (org.EcK12.eg.db), E. coli strain Sakai (org.EcSakai.eg.db), zebrafish (org.Dr.eg.db) and Drosophila (org.Dm.eg.db) were installed with the required R packages. If your species is not in the above list, contact the eResearch team.

## USER INPUT
# Input your species genome below (select from the above list)
my_genome <- org.Mm.eg.db

# Pull out just the Entrez/Ensembl gene IDs
gene_ids <- row.names(res)

# You can see the list of annotations you can apply to your data by:
keytypes(my_genome)

# Annotate gene symbol and description to your DE gene IDs (you can add other keytypes from the above 'keytypes(my_genome)' list, if you choose)
cols <- c("SYMBOL", "GENENAME")

## USER INPUT
# Provide the gene ID type for your DE data
# If you have Ensemble IDs, enter "ENSEMBL", if you have Entrez IDs, enter "ENTREZID", if you have gene symbols, enter "SYMBOL"
idtype <- "ENSEMBL"

# Map the Entrez/Ensembl gene IDs to gene symbol and description
map <- AnnotationDbi::select(my_genome, keys=gene_ids, columns=cols, keytype=idtype)

# Combine the annotation data with the DE data
# Since we're matching Ensembl -> Entrez there aren't a 1:1 mapping, so need to merge rather than cbind
annot <- merge(x = res, y = map, by.x = 0, by.y = idtype, all = F)
# There isn't always a 1:1 mapping between gene identifiers, so we also need to remove duplicates
annot_nodups <- annot[!duplicated(annot$Row.names), ]
# Reorder by adjusted p
annot_nodups <- annot_nodups[order(annot_nodups$padj), ]

# Add normalised counts to the output table. This is so you can later plot expression trends for individual genes in R, Excel, etc.
# Need to normalise the counts first, using the size factors calculated by DESeq2 (in the 'deseq' object)
expdata_norm <- as.matrix(expdata) %*% diag(deseq$sizeFactor)
colnames(expdata_norm) <- colnames(expdata)
annot_counts <- merge(x = annot_nodups, y = expdata_norm, by.x = "Row.names", by.y = 0, all = TRUE)

# Export the full annotated dataset as a csv (for journal submission)
# Creates a 'Tables' subdirectory where all tables will be output
dir.create("Tables", showWarnings = FALSE)
write.csv(annot_counts, file=paste0("./Tables/all_genes_", paste(degroups, collapse = "_Vs_"), ".csv"), row.names = FALSE)

# Pull out just significant genes (change from 0.05 to 0.01 if needed)
DE_genes <- subset(annot_counts, padj < 0.05, select=colnames(annot_counts))
# NOTE: add lfc cutoffs if needed. E.g., log2FoldChange > 1 and < log2FoldChange -1 cutoff
# DE_genes <- subset(DE_genes, log2FoldChange > 1 | log2FoldChange < -1, select=colnames(DE_genes))

# Export the sig DE genes as a table (to be used in downstream analysis)
write.csv(DE_genes, file=paste0("./Tables/DE_genes_", paste(degroups, collapse = "_Vs_"), ".csv"), row.names = FALSE)


#### 6b. Volcano plot ####

p <- EnhancedVolcano(annot_nodups, lab = annot_nodups$SYMBOL, selectLab = annot_nodups$SYMBOL[1:20], drawConnectors = TRUE, title = NULL, subtitle = NULL, x = 'log2FoldChange', y = 'pvalue')

p

# NOTE: the above plot shows labels for the top significantly DE (i.e. by lowest adjusted p value) genes.

# Output as publication quality (300dpi) tiff and pdf.
# Create a (300dpi) tiff
ggsave(file = paste0("./figures/volcano_", paste(degroups, collapse = "_Vs_"), ".tiff"), dpi = 300, compression = "lzw", device = "tiff", width = 10, height = 8, plot = p)
# Create a pdf
ggsave(file = paste0("./figures/volcano_", paste(degroups, collapse = "_Vs_"), ".pdf"), device = "pdf", width = 10, height = 8, plot = p)


#### 6c. DE genes heatmaps and dendrograms ####

# Make the row names gene symbols.
DE_genes <- na.omit(DE_genes)
row.names(DE_genes) <- make.unique(DE_genes$SYMBOL)

# sort by p-value
DE_genes <- DE_genes[order(DE_genes$padj), ]

# Pull out normalised counts only
siggc <- DE_genes[colnames(DE_genes) %in% colnames(expdata)]

# Scale and center each row. This is important to visualise relative differences between groups and not have row-wise colouration dominated by high or low gene expression.
xts <- scale(t(siggc))
xtst <- t(xts)

# Define annotation column
annot_columns <- data.frame(meta$group[meta$group %in% degroups])
# Make the row names the sample IDs
row.names(annot_columns) <- meta$sample_ID[meta$group %in% degroups]
colnames(annot_columns) <- "Treatment groups"
# Need to factorise it
annot_columns[[1]] <- factor(annot_columns[[1]])

# Generate dendrogram and heatmap for ALL DE genes
pheatmap(xtst, color=colorRampPalette(c("#D55E00", "white", "#0072B2"))(100), annotation_col=annot_columns, annotation_names_col = F, fontsize_col = 12, fontsize_row = 7, labels_row = row.names(siggc), show_rownames = F, filename = paste0("./figures/All_DEG_Heatmap_", paste(plotgroups, collapse = "_Vs_"), ".tiff"))

# Generate dendrogram and heatmap for top 20 DE genes
pheatmap(xtst[1:20,], color=colorRampPalette(c("#D55E00", "white", "#0072B2"))(100), annotation_col=annot_columns, annotation_names_col = F, fontsize_col = 14, fontsize_row = 12, labels_row = row.names(siggc), show_rownames = T, filename = paste0("./figures/Top_DEG_Heatmap_", paste(plotgroups, collapse = "_Vs_"), ".tiff"))
# NOTE: you can plot more than 20 top genes by adjusting 'xtst[1:20,]'. If you wanted to plot the top 50 genes you'd change this to 'xtst[1:50,]'




