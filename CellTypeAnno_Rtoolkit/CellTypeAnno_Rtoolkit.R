print("=========Task Started=========")

suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleR)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

args <- commandArgs(trailingOnly=TRUE)
ref_dir <- args[1]
query_dir <- args[2]
out_dir <- args[3]

print("=========Build Reference Dataset=========")
ref_obj <- readRDS(ref_dir)
ref_obj <- NormalizeData(ref_obj)
ref_obj <- FindVariableFeatures(ref_obj, selection.method = "vst", nfeatures = 2000)

ref_counts <- as.matrix(ref_obj@assays$RNA@data[VariableFeatures(ref_obj),])
label <- ref_obj$cell_type
ref_se <- SummarizedExperiment(assays=list(counts = ref_counts, logcounts = ref_counts))
ref_se@colData$label.main <- label
print("Finished")

print("=========Annotation Started=========")
query_obj <- readRDS(query_dir)
query_obj <- NormalizeData(query_obj)
query_obj <- FindVariableFeatures(query_obj, selection.method = "vst", nfeatures = 2000)

query_counts <- as.matrix(query_obj@assays$RNA@data[VariableFeatures(query_obj),])
query_sce <- SingleCellExperiment(assays=list(counts=query_counts, logcounts = query_counts))

query_pred <- SingleR(test = query_sce, ref = ref_se, labels = label)
query_result <- data.frame(main_type=query_pred$labels)
rownames(query_result) <- colnames(query_sce)
print("Finished")

print("=========Writing Output File=========")
write.csv(query_result, paste0(query_dir,".csv"))

