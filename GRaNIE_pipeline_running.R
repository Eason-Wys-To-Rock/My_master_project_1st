##--------------------GRaNIE pipeline running--------------------
#Set working directory, file name, and file paths
work_dir = "/exports/eddie/scratch/s2858981/data_GRN_test/GRaNIE_formal_test/GRN_control/"
output_filename = "/GRN_atac_all_rna_1w_sam_control_clustered.rds"
rna_csv = "/exports/eddie/scratch/s2858981/data/output/cluster_and_pseudobulk/sam_clustered/sam_rna_df_agg_control.csv"
atac_csv = "/exports/eddie/scratch/s2858981/data/output/cluster_and_pseudobulk/sam_clustered/sam_atac_df_agg_control.csv"
metadata_csv = "/exports/eddie/scratch/s2858981/data/output/cluster_and_pseudobulk/sam_clustered/sam_metadata_df_agg_control.csv"

#Library import
print("----------Packages recruiting start----------") # start recruiting
library(readr)
library(GRaNIE)
library(org.Mm.eg.db) # This is for mouse
print("----------Packages recruiting finish----------") # end recruiting

#import values
idColumn_peaks = "peakID"
idColumn_RNA = "ENSEMBL"
genomeAssembly = "mm10"
dir_output = "."

##--------------------GRaNIE start running--------------------
print("Now start running pipeline for Control Group")
#Import dataset
setwd(work_dir)
GRN_file_outputRDS = paste0(dir_output, output_filename)
print("Basic parameters are prepared. Including idColumn_peaks, idColumn_RNA, genomeAssembly, and dir_output.")

# Introduce RNA
rna_df <- read.csv(rna_csv)
# Introduce ATAC
atac_df <- read.csv(atac_csv)
#Introduce metadata
metadata_df <- read.csv(metadata_csv)
rna_df$X <- NULL
atac_df$X <- NULL
metadata_df$X <- NULL
print("Data Introduce Done.")

# Save the name of the respective ID columns. Note: for rna, in GRaNIE it was "ENSEMBL" format, but not in our data.
GRN = initializeGRN(outputFolder = dir_output, genomeAssembly = genomeAssembly)
message("GRN is initialized")

# Add Data
#This may take a bit of time...
GRN = addData(GRN, counts_peaks = atac_df, normalization_peaks = "none",
              idColumn_peaks = idColumn_peaks, counts_rna = rna_df, normalization_rna = "none",
              idColumn_RNA = idColumn_RNA, sampleMetadata = metadata_df, forceRerun = TRUE)
message("Data is added to GRN")

#delete values to save memory
rm(rna_df, atac_df, metadata_df)

#add TF
if (dir.exists("PWMScan_HOCOMOCOv12")) {
  
  # if TRUE, then continue
  message("TFData directory found. Continue running...")
  
} else {
  
  # if FALSE, stop running and return ERROR
  stop("TFData is not in the working directory")
}
motifFolder = tools::file_path_as_absolute("PWMScan_HOCOMOCOv12")

GRN = addTFBS(GRN, motifFolder = motifFolder, TFs = "all", filesTFBSPattern = "_TFBS",
              fileEnding = ".bed.gz", forceRerun = TRUE)
message("TFBS is added")

GRN = overlapPeaksAndTFBS(GRN, nCores = 1, forceRerun = TRUE)
message("TFBS is overlapped to peaks")

#Save Spot 1
message("Save Spot 1 start...")
#fix saveRDS error after initialize GRN.
setMethod("containsOutOfMemoryData", "GRN", function(x) FALSE)
#saveRDS
saveRDS(GRN, GRN_file_outputRDS)
message("Save Spot 1 is done")
print("Save Spot 1 is done")


#--------------------Filter-Data-Step-is-excluded--------------------
message("Filter-Data-Step-is-excluded")
#--------------------Add-Connections--------------------
message("Add-Connections-Start")
#Add TF-Enhancer Connections
GRN = addConnections_TF_peak(GRN, plotDiagnosticPlots = FALSE, connectionTypes = c("expression"),
                             corMethod = "spearman", forceRerun = TRUE)

#Add Enhancer-gene Connections
GRN = addConnections_peak_gene(GRN, corMethod = "spearman", promoterRange = 10000,
                               TADs = NULL, nCores = 1, plotDiagnosticPlots = FALSE, plotGeneTypes = list(c("all")),
                               forceRerun = TRUE)
#Combine TF-enhancer and enhancer-gene connections and filter
GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = 0.2, peak_gene.fdr.threshold = 0.2,
                               peak_gene.fdr.method = "BH", gene.types = c("protein_coding", "lincRNA"), allowMissingTFs = FALSE,
                               allowMissingGenes = FALSE, forceRerun = TRUE)
print("TF-peak and peak-gene fdr thersholds are set to 0.2")

#(optional But I did it) Add TF-gene correlations
GRN = add_TF_gene_correlation(GRN, corMethod = "spearman", nCores = 1, forceRerun = TRUE)

#Save Spot 2
message("Save Spot 2 start...")
#fix saveRDS error after initialize GRN.
saveRDS(GRN, GRN_file_outputRDS)
message("Save Spot 2 is done")
print("Save Spot 2 is done")

#--------------------Retrieve-filtered-connections--------------------
GRN_connections.all = getGRNConnections(GRN, type = "all.filtered", include_TF_gene_correlations = TRUE,
                                        include_geneMetadata = TRUE)

print("GRN_connections.all")
GRN_connections.all

#--------------------Construct-GRN--------------------
GRN = build_eGRN_graph(GRN, forceRerun = TRUE)

#Save the Final GRN.
message("Save Final start...")
saveRDS(GRN, GRN_file_outputRDS)
message("Save Final is done")
print("Save Final is done")

#--------------------Start-Plotting--------------------
#QC2: Diagnostic plots for TF-enhancer connections
#GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("real"), plotAsPDF = FALSE, pages = c(1))
#GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("real"), plotAsPDF = FALSE, pages = c(42))
#GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("background"), plotAsPDF = FALSE, pages = c(42))

#QC3:
#GRN = plotDiagnosticPlots_peakGene(GRN, gene.types = list(c("protein_coding", "lincRNA")), plotAsPDF = FALSE, pages = 1)

#Generate a Connection Summary for Filtered Connections
GRN = generateStatsSummary(GRN, TF_peak.fdr = c(0.05, 0.1, 0.2), TF_peak.connectionTypes = "all",
                           peak_gene.fdr = c(0.1, 0.2), peak_gene.r_range = c(0, 1),
                           allowMissingGenes = c(FALSE, TRUE), allowMissingTFs = c(FALSE),
                           gene.types = c("protein_coding", "lincRNA"),
                           forceRerun = TRUE)
GRN = plot_stats_connectionSummary(GRN, type = "heatmap", plotAsPDF = TRUE, pages = 3)

#Visualize the GRN
GRN = visualizeGRN(GRN, plotAsPDF = TRUE, maxEdgesToPlot = 1000)

#Generate Network and enrichment analyses for filtered connections
GRN = performAllNetworkAnalyses(GRN, ontology = c("GO_BP"), outputFolder = ".", forceRerun = TRUE)

GRN
#Save the Final GRN again.
saveRDS(GRN, GRN_file_outputRDS)
rm(GRN, GRN_file_outputRDS, GRN_connections.all, motifFolder)
gc()
message("GRaNIE analysis is done.")
print("GRaNIE analysis is done.")
