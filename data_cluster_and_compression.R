### This code aggregate sam's atac&rna data, then cluster & do the pseudobulk expression for GRaNIE's further analysis
#This data doesn't have file checking steps. Be careful if you want to use it!
message("This data doesn't have file checking steps. Be careful if you want to use it!")

# set working directory and file read path
setwd("/exports/eddie/scratch/s2858981/data/output/cluster_and_pseudobulk/")
seurat_rna_rds <- "/exports/eddie/scratch/s2858981/data/input/scenic_scRNA_data/seurat-mm10-rna.rds"
seurat_atac_rds <- "/exports/eddie/scratch/s2858981/data/input/scenic_scATAC/seurat-mm10-atac.rds"

#Setting the genotype: Control or Pax6-cKO. This will be used for following data integration and compression.
subset_genotype = "Control"

###----------file save pathway----------
multiome.file <- "./sam_reduction_done_control.rds"
metadata.file <- "./sam_clustered/sam_metadata_df_agg_control.csv"
rna.file <- "./sam_clustered/sam_rna_df_agg_control.csv"
atac.file <- "./sam_clustered/sam_atac_df_agg_control.csv"

###--------------------Package import--------------------
library(Seurat)
library(Signac)
library(biomaRt)
library(dplyr)
library(tibble)

###Function setting
##Function 1
#import add column function

add_id_column <- function(df, id_col_name = "ID") {
  # check
  if (!is.data.frame(df)) stop("Input must be a data.frame")
  if (is.null(rownames(df))) stop("Input data.frame has no rownames")
  
  # Put rownames into the 1st column
  df_new <- cbind(setNames(list(rownames(df)), id_col_name), df)
  
  # delete the rownames
  rownames(df_new) <- NULL
  
  return(df_new)
}

##function 2
#convert the ATAC ID in dense matrix. This only changes rowname. Note: This only works in Kai's dataset for GRaNIE analysis.
# Convert ATAC peak names in the rownames of a dense matrix
fix_peak_rownames_dense <- function(mat) {
  rn <- rownames(mat)
  
  # split the name into 3 pieces by "-" in the text
  parts <- do.call(rbind, strsplit(rn, "-", fixed = TRUE))
  
  # create new rownames
  new_rn <- paste0("chr", parts[,1], ":", parts[,2], "-", parts[,3])
  
  # write back rownames，without changing other dataframe information.
  rownames(mat) <- new_rn
  
  return(mat)
}
### End of function definitions




###----------subsetting by genotype if needed----------
#multiome <- subset(
#  x = multiome,
#  subset = genotype == "Control"
#)

#import Seurat data
seurat_rna <- readRDS(seurat_rna_rds)
seurat_atac <- readRDS(seurat_atac_rds)

#subsetting data by genotype if needed
seurat_rna_sub <- subset(
  x = seurat_rna,
  subset = genotype == subset_genotype
)
seurat_atac_sub <- subset(
  x = seurat_atac,
  subset = genotype == subset_genotype
)

#Find Common cells
common_cells  <- colnames(multiome)
#check the length
length(common_cells)
#Subsetting common cell barcodes and creat assay object suitable for seurat analysis
seurat_rna_sub <- subset(seurat_rna, cells = common_cells)
seurat_atac_sub <- subset(seurat_atac, cells = common_cells)
#Find Variable/Top features
seurat_rna_sub <- FindVariableFeatures(seurat_rna_sub, selection.method = "vst", nfeatures = 10000)
multi_variablefeatures <- VariableFeatures(seurat_rna_sub)
#(Optional step at the following. You can run it if the atac information is too much/)
#seurat_atac_sub <- FindTopFeatures(seurat_atac_sub, min.cutoff = "q50")
#multi_topfeatures <- rownames(seurat_atac_sub@assays$originalexp@counts)[seurat_atac_sub@assays$originalexp@meta.features$percentile >= 0.5]

#transfer the data into original project.
multiome@assays$RNA@meta.features <- seurat_rna_sub@assays$originalexp@meta.features
multiome@assays$ATAC@meta.features <- seurat_atac_sub@assays$originalexp@meta.features

#remove items
rm(seurat_rna, seurat_rna_sub, seurat_atac, seurat_atac_sub)
gc()


# Perform standard analysis of each modality independently RNA analysis
multiome <- FindVariableFeatures(multiome)
multiome <- ScaleData(multiome)
multiome <- RunPCA(multiome)
multiome <- RunUMAP(multiome, dims = 1:50, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
message("RNA reduction is done")
# ATAC analysis add gene annotation information
multiome <- FindTopFeatures(multiome, min.cutoff = "q0")
multiome <- RunSVD(multiome)
multiome <- RunUMAP(multiome, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
message("ATAC reduction is done")
message("both reduction is done!!!")
print("RNA and ATAC reduction is done")

multiome <- FindMultiModalNeighbors(multiome, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
multiome <- RunUMAP(multiome, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
multiome <- FindClusters(multiome, graph.name = "wsnn", algorithm = 4, resolution = 2, verbose = FALSE)
#save file, rename if required.
saveRDS(multiome, multiome.file)
print("clustering finished. Now start aggregation...")


multiome_metadata <- multiome@meta.data %>%
  dplyr::select(
    cluster = "wsnn_res.2"
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::arrange(as.numeric(as.character(cluster))) %>%
  dplyr::rename(
    cluster_number = cluster
  )
#Further groups are started with numbers, and AggregateExpression() appends "g" before the numbers. So paste g before numbers in metadata manually.
multiome_metadata[[1]] <- paste0("g", multiome_metadata[[1]])
#rename the metadata[[1]] into "sample_id"
colnames(multiome_metadata)[1] <- "sample_id"

#Metadata output, rename if required.
write.csv(multiome_metadata, metadata.file)

##Export the clustered datasets
agg_multiome_whole <- AggregateExpression(
  multiome,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "wsnn_res.2",
  add.ident = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE,
)

#extract RNA (Note: use the correct dataframe!)
rna_df_agg <- as.data.frame(agg_multiome_whole$RNA)
#gene filter by variable features
rna_df_agg <- rna_df_agg[multi_variablefeatures, ]

#create column with ENSEMBL volcabulary
#> packageVersion("biomaRt")
#[1] ‘2.66.0’
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

gene2ens <- getBM(
  attributes = c(
    "external_gene_name",
    "ensembl_gene_id"
  ),
  filters = "external_gene_name",
  values = rownames(rna_df_agg),
  mart = mart
)

#format changing: requirements:
#1. change colname from ensembl_gene_id to ENSEMBL so that it can be recognized in GRaNIE.
#2. delete the column "X" which is unnecessary.
colnames(gene2ens)[2] <- "ENSEMBL"

gene2ens <- gene2ens %>%
  distinct(external_gene_name, .keep_all = TRUE) %>%
  column_to_rownames(var = "external_gene_name")

ensembl_vec <- gene2ens[rownames(rna_df_agg), 1]

rna_df_agg$ENSEMBL <- ensembl_vec

rna_df_agg <- rna_df_agg[, c("ENSEMBL", setdiff(colnames(rna_df_agg), "ENSEMBL"))]

rna_df_agg <- rna_df_agg[!is.na(rna_df_agg$ENSEMBL), ]

rownames(rna_df_agg) <- NULL

write.csv(rna_df_agg, rna.file)

## ATAC (Note: use the correct dataframe!)
atac_df_agg <- as.data.frame(agg_multiome_whole$ATAC)
#ATAC filter by top features
#atac_df_agg <- atac_df_agg[rownames(atac_df_agg) %in% multi_topfeatures, ]

atac_df_agg <- add_id_column(atac_df_agg, id_col_name = "peakID")

write.csv(atac_df_agg, atac.file)

rm(atac_df_agg, rna_df_agg, agg_multiome_whole, multiome_metadata, multi_variablefeatures, multi_topfeatures, multiome, mart, gene2ens, ensembl_vec)
gc()

print("data processing is done. Now for GRaNIE pipeline analysis")
message("data processing is done. Now for GRaNIE pipeline analysis")