library("CHETAH")
library("scCATCH")
library("scmap")
library("SingleR")
library("Seurat")
library("stringr")
library("scater")
library('SCINA')
library('preprocessCore')

# design for batch runnning
query_dir_name_list <- c(

"/stor/public/hcad/brain_cortex_Gaublomme2019_part1/brain_cortex_Gaublomme2019_part1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_cortex_Gaublomme2019_part2/brain_cortex_Gaublomme2019_part2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_MTG_AllenBrainAtlas/brain_MTG_AllenBrainAtlas.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_VisualCortexFrontalCortexCerebellum_Lake2017/brain_VisualCortexFrontalCortexCerebellum_Lake2017.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_rett_William2018/brain_rett_William2018.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_hippocampus_Zhong2020/brain_hippocampus_Zhong2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_ACC_AllenBrainAtlas/brain_ACC_AllenBrainAtlas.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_V1_AllenBrainAtlas/brain_V1_AllenBrainAtlas.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_PFC_Zhong2018/brain_PFC_Zhong2018.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/eye_retina_Lukowski2019/eye_retina_Lukowski2019.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_ParietalTemporalFrontalFrontoparietal_Venteicher2017/brain_ParietalTemporalFrontalFrontoparietal_Venteicher2017.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_TemporalLobe_Guo2020/brain_TemporalLobe_Guo2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_AdultCerebellum1_Guo2020/brain_AdultCerebellum1_Guo2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_FetalBrain3_Guo2020/brain_FetalBrain3_Guo2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_FetalBrain4_Guo2020/brain_FetalBrain4_Guo2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_FetalBrain5_Guo2020/brain_FetalBrain5_Guo2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/brain_FetalBrain6_Guo2020/brain_FetalBrain6_Guo2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Transverse-Colon_Transverse-Colon_HCLAdultTransverse-Colon2/Transverse-Colon_Transverse-Colon_HCLAdultTransverse-Colon2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Sigmoid-Colon_Sigmoid-Colon_HCLAdultSigmoid-Colon1/Sigmoid-Colon_Sigmoid-Colon_HCLAdultSigmoid-Colon1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine1/Intestine_Intestine_HCLFetalIntestine1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine2/Intestine_Intestine_HCLFetalIntestine2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine3/Intestine_Intestine_HCLFetalIntestine3.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine4/Intestine_Intestine_HCLFetalIntestine4.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine5/Intestine_Intestine_HCLFetalIntestine5.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Rectum_Rectum_HCLAdultRectum1/Rectum_Rectum_HCLAdultRectum1.seuratobj.dbupload_v1.rds",

"/stor/public/hcad/lung_lung_Madissoon2019/lung_lung_Madissoon2019.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/lung_lung_carcinomas2019/lung_lung_carcinomas2019.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/heart_heart_asp2019cellpress/heart_heart_asp2019cellpress.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/heart_heart_HCLAdultHeart1/heart_heart_HCLAdultHeart1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/heart_heart_HCLAdultHeart2/heart_heart_HCLAdultHeart2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/heart_heart_HCLFetalHeart1/heart_heart_HCLFetalHeart1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/heart_heart_HCLFetalHeart2/heart_heart_HCLFetalHeart2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Heart_Heart_Wang2020/Heart_Heart_Wang2020.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Heart_Heart_Wang2020_2/Heart_Heart_Wang2020_2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/lung_lung_HCLAdultLung1/lung_lung_HCLAdultLung1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/lung_lung_HCLAdultLung2/lung_lung_HCLAdultLung2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/lung_lung_HCLAdultLung3/lung_lung_HCLAdultLung3.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/lung_lung_HCLFetalLung1/lung_lung_HCLFetalLung1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/lung_lung_HCLFetalLung2/lung_lung_HCLFetalLung2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Blood_Blood_HCLCord-Blood1/Blood_Blood_HCLCord-Blood1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Blood_Blood_HCLCord-Blood2/Blood_Blood_HCLCord-Blood2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood1/Blood_Blood_HCLAdult-Peripheral-Blood1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood2/Blood_Blood_HCLAdult-Peripheral-Blood2.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood3/Blood_Blood_HCLAdult-Peripheral-Blood3.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood4/Blood_Blood_HCLAdult-Peripheral-Blood4.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow1/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/Spleen_Spleen_HCLAdult-Spleen1/Spleen_Spleen_HCLAdult-Spleen1.seuratobj.dbupload_v1.rds",
"/stor/public/hcad/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow2/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow2.seuratobj.dbupload_v1.rds"
)


# default query dir "/stor/public/hcad/", change it to your dir
# query_dir_list <- paste("/stor/public/hcad/", query_dir_name_list, "/", query_dir_name_list, ".seuratobj.rds", sep = "")
query_dir_list <- query_dir_name_list

From_hcadTissue_To_sccatchTissue <- function(hcadTissue){
  # transform to lowercase
  hcadTissue = tolower(hcadTissue)
  
  if(hcadTissue == "prostate"){
    return("Prostate")
  }else if(hcadTissue == "mammary gland"){
    return("Mammary gland")
  }else if(hcadTissue == "breast"){
    return("Breast")
  }else if(hcadTissue == "liver"){
    return("Liver")
  }else if(hcadTissue == "gallbladder"){
    return("Gall bladder")
  }else if(hcadTissue == "pancreas"){
    return("Pancreas")
  }else if(hcadTissue == "adrenal gland"){
    return("Adrenal gland")
  }else if(hcadTissue == "thyroid gland"){
    return("Thyroid")
  }else if(hcadTissue == "heart"){
    return("Heart")
  }else if(hcadTissue == "vein"){
    return("Blood vessel")
  }else if(hcadTissue == "artery"){
    return("Blood vessel")
  }else if(hcadTissue == "lung"){
    return("Lung")
  }else if(hcadTissue == "nasal airway"){
    return("Airway epithelium")
  }else if(hcadTissue == "bronchus"){
    return("Bronchiole")
  }else if(hcadTissue == "trachea"){
    return("Trachea")
  }else if(hcadTissue == "nose"){
    return("Nasal concha")
  }else if(hcadTissue == "larynx"){
    return("Larynx")
  }else if(hcadTissue == "pleura"){
    return("Pleura")
  }else if(hcadTissue == "oesophagus"){
    return("Esophagus")
  }else if(hcadTissue == "intestine"){
    return("Intestine")
  }else if(hcadTissue == "stomach"){
    return("Stomach")
  }else if(hcadTissue == "kidney"){
    return("Kidney")
  }else if(hcadTissue == "bladder"){
    return("Bladder")
  }else if(hcadTissue == "ureter"){
    return("Bladder")
  }else if(hcadTissue == "urethra"){
    return("Bladder")
  }else if(hcadTissue == "adipose"){
    return("Adipose tissue")
  }else if(hcadTissue == "muscle"){
    return("Muscle")
  }else if(hcadTissue == "cartilage"){
    return("Cartilage")
  }else if(hcadTissue == "osseous"){
    return("Bone")
  }else if(hcadTissue == "rib"){
    return("Bone")
  }else if(hcadTissue == "skin"){
    return("Skin")
  }else if(hcadTissue == "blood"){
    return("Blood")
  }else if(hcadTissue == "bone marrow"){
    return("Bone marrow")
  }else if(hcadTissue == "tonsil"){
    return("Tonsil")
  }else if(hcadTissue == "immune cell"){
    return("Lymphoid tissue")
  }else if(hcadTissue == "lymph"){
    return("Lymph")
  }else if(hcadTissue == "spleen"){
    return("Spleen")
  }else if(hcadTissue == "thymus"){
    return("Thymus")
  }else if(hcadTissue == "pns"){
    return("Sympathetic ganglion")
  }else if(hcadTissue == "neuron/glia"){
    return("Brain")
  }else if(hcadTissue == "brain"){
    return("Brain")
  }else if(hcadTissue == "spinal cord"){
    return("Spinal cord")
  }else if(hcadTissue == "retina"){
    return("Retina")
  }else if(hcadTissue == "eye"){
    return("Eye")
  }else if(hcadTissue == "uterus"){
    return("Uterus")
  }else if(hcadTissue == "cervix"){
    return("Uterus")
  }else if(hcadTissue == "fallopian Tube"){
    return("Ovary")
  }else if(hcadTissue == "testis"){
    return("Testis")
  }else if(hcadTissue == "ovary"){
    return("Ovary")
  }else if(hcadTissue == "embryo"){
    return("Embryo")
  }else if(hcadTissue == "placenta"){
    return("Placenta")
  }else{
    return("UNKNOWN")
  }
}

# main structure
for(dataset in 1:length(query_dir_name_list)){
  print(paste("=======Querying dataset ", query_dir_name_list[dataset], " ========"))
  # load query data
  query_obj <- readRDS(query_dir_list[dataset])
  
  # The dataframe to save the annotation results
  annotation_results <- data.frame(Tissue = as.character(query_obj$organ))
  rownames(annotation_results) <- rownames(query_obj@meta.data)
  
  query_obj <- FindVariableFeatures(query_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(query_obj)
  query_obj <- ScaleData(query_obj, features = all.genes)
  query_obj <- RunPCA(query_obj, features = VariableFeatures(object = query_obj))
  
  # scCATCH
  print("====Running scCATCH====")
  query_obj_sccatch <- query_obj
  query_obj_sccatch <- FindNeighbors(query_obj_sccatch)
  query_obj_sccatch <- FindClusters(query_obj_sccatch, resolution = 0.5)
  
  clu_markers <- findmarkergenes(object = query_obj_sccatch,
                                 species = 'Human',
                                 cluster = 'All',
                                 match_CellMatch = FALSE,
                                 cancer = NULL,
                                 tissue = NULL,
                                 cell_min_pct = 0.25,
                                 logfc = 0.25,
                                 pvalue = 0.05)
  query_tissue <- From_hcadTissue_To_sccatchTissue(as.character(query_obj_sccatch$organ[1]))
  
  if(query_tissue != "UNKNOWN"){
    clu_ann <- scCATCH(object = clu_markers$clu_markers,
                       species = 'Human',
                       cancer = NULL,
                       tissue = query_tissue)
    cell_ident <- as.numeric(query_obj_sccatch@active.ident)
    query_type <- as.numeric(query_obj_sccatch@active.ident)
    for (i in c(1:length(query_type))){
      query_type[i] <- clu_ann$cell_type[cell_ident[i]]
    }
    query_result <- data.frame(cell_type=query_type)
    rownames(query_result) <- colnames(query_obj_sccatch@assays$RNA@counts)
    
    ## default result dir, your local folder
    # save_dir <- paste(query_dir_name_list[dataset], "_predictions_scCATCH.csv", sep = "")
    # write.csv(query_result, file = save_dir)
    
    ## clear environment
    rm(query_obj_sccatch, clu_markers, clu_ann)
    
    ## save in annotation_results
    annotation_results$scCATCH_annotation <- as.character(query_result[,1])
  }
  
  ## load the shared ref for scmap, here we use HPCA for all query data
  ## This treatment is specially designed for HPCA, in later versions may change
  ref_obj <- readRDS("/stor/public/share/HPCA/HPCA.seuratobj.rds")
  counts.mat <- as.matrix(ref_obj@assays$RNA@counts)
  sample.size <- Matrix::colSums(counts.mat)
  for(i in 1:745){
    counts.mat[,i] = counts.mat[,i]/sample.size[i]*10^6
  }
  logcounts.mat <- log2(counts.mat + 1)
  label = ref_obj$hpca.main.label
  
  ## chetah
  print("====Running CHETAH====")
  ref_sce <- SingleCellExperiment(assays = list(counts = counts.mat),colData = DataFrame(celltypes = label))
  query_sce <- SingleCellExperiment(assays=list(counts = query_obj@assays$RNA@counts))
  input <- CHETAHclassifier(input = query_sce, ref_cells = ref_sce)
  query_result <- input$celltype_CHETAH
  
  ## default result dir, your local folder
  # save_dir <- paste(query_dir_name_list[dataset], "_predictions_CHETAH.csv", sep = "")
  # write.csv(query_result, file = save_dir)
  
  ## save in annotation_results
  annotation_results$CHETAH_annotation <- as.character(query_result)
  
  ## clear environment
  rm(ref_sce, query_sce, input)
  
  ## scmap
  print("====Running scmap====")
  ref_sce <- SingleCellExperiment(assays=list(counts = counts.mat, logcounts = logcounts.mat))
  ref_sce@colData$cell_type1 <- label
  rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
  ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]
  ### selectFeatures may return error, I recommend setFeatures
  var_genes <- VariableFeatures(query_obj)
  var_genes <- str_replace(var_genes, "-", "_")
  ref_sce <- setFeatures(ref_sce, var_genes)
  ref_sce <- indexCluster(ref_sce)
  
  query_sce <- SingleCellExperiment(assays=list(counts=as.matrix(query_obj@assays$RNA@counts)))
  logcounts(query_sce) <- log2(counts(query_sce) + 1)
  rowData(query_sce)$feature_symbol <- rownames(query_sce)
  
  scmapCluster_results <- scmapCluster(projection = query_sce, index_list = list(ref = metadata(ref_sce)$scmap_cluster_index))
  query_result <- data.frame(main_type=scmapCluster_results$scmap_cluster_labs)
  
  ## default result dir, your local folder
  # save_dir <- paste(query_dir_name_list[dataset], "_predictions_scmap.csv", sep = "")
  # write.csv(query_result, file = save_dir)
  
  ## save in annotation_results
  annotation_results$scmap_annotation <- as.character(query_result[,1])
  
  ## clear environment
  rm(ref_sce, query_sce, scmapCluster_results)
  
  # singler
  print("====Running SingleR====")
  ref_se <- SummarizedExperiment(assays=list(counts = counts.mat, logcounts = logcounts.mat))
  ref_se@colData$label.main <- label
  query_sce <- SingleCellExperiment(assays=list(counts=query_obj@assays$RNA@counts))
  query_sce <- logNormCounts(query_sce)
  
  query_pred <- SingleR(test = query_sce, ref = ref_se, labels = label)
  query_result <- data.frame(main_type=query_pred$labels)
  rownames(query_result) <- colnames(query_sce)
  
  ## default result dir, your local folder
  # save_dir <- paste(query_dir_name_list[dataset], "_predictions_SingleR.csv", sep = "")
  # write.csv(query_result, file = save_dir)
  
  ## save in annotation_results
  annotation_results$SingleR_annotation <- as.character(query_result[,1])
  
  ## clear environment
  rm(ref_se, query_sce, query_pred)
  
  # scina
  print("====Running scINA====")
  ## load signatures of cell types (prior knowlege)
  ## prior knowledge: signatures (A list contains multiple signature vectors)
  ## using eTME signatures
  load(system.file('extdata','example_signatures.RData', package = "SCINA"))
  ## using own signatures files:
  ## signatures=preprocess.signatures('your/path/to/example_signatures.csv')
  
  ## input data: expression matrix 
  ## From .rds
  exp <- query_obj@assays$RNA@scale.data
  
  ## running scina
  results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
                  convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
  query_result <- data.frame(results$cell_labels)
  
  ## default result dir, your local folder
  # save_dir <- paste(query_dir_name_list[dataset], "_predictions_SCINA.csv", sep = "")
  # write.csv(query_result, file = save_dir)
  
  ## save in annotation_results
  annotation_results$SCINA_annotation <- as.character(query_result[,1])
  
  ## clear environment
  rm(exp, results)
  
  # save the merged annotation result, default in your local folder
  save_dir <- paste(unlist(strsplit(query_dir_name_list[dataset],'/'))[5], "_predictions.rds", sep = "")
  saveRDS(annotation_results, file = save_dir)
}


# # Total list
# query_dir_name_list <- c(
#   "/stor/public/hcad/brain_LGN_AllenBrainAtlas/brain_LGN_AllenBrainAtlas.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/eye_retina_Menon2019_10x/eye_retina_Menon2019_10x.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/eye_retina_Menon2019_seq-well/eye_retina_Menon2019_seq-well.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/retina_macula_Voigt2019/retina_macula_Voigt2019.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_cortex_Gaublomme2019_part1/brain_cortex_Gaublomme2019_part1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_cortex_Gaublomme2019_part2/brain_cortex_Gaublomme2019_part2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_MTG_AllenBrainAtlas/brain_MTG_AllenBrainAtlas.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_VisualCortexFrontalCortexCerebellum_Lake2017/brain_VisualCortexFrontalCortexCerebellum_Lake2017.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_rett_William2018/brain_rett_William2018.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_hippocampus_Zhong2020/brain_hippocampus_Zhong2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_ACC_AllenBrainAtlas/brain_ACC_AllenBrainAtlas.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_V1_AllenBrainAtlas/brain_V1_AllenBrainAtlas.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_PFC_Zhong2018/brain_PFC_Zhong2018.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/eye_retina_Lukowski2019/eye_retina_Lukowski2019.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_ParietalTemporalFrontalFrontoparietal_Venteicher2017/brain_ParietalTemporalFrontalFrontoparietal_Venteicher2017.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_TemporalLobe_Guo2020/brain_TemporalLobe_Guo2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_AdultCerebellum1_Guo2020/brain_AdultCerebellum1_Guo2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_FetalBrain3_Guo2020/brain_FetalBrain3_Guo2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_FetalBrain4_Guo2020/brain_FetalBrain4_Guo2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_FetalBrain5_Guo2020/brain_FetalBrain5_Guo2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/brain_FetalBrain6_Guo2020/brain_FetalBrain6_Guo2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Transverse-Colon_Transverse-Colon_HCLAdultTransverse-Colon2/Transverse-Colon_Transverse-Colon_HCLAdultTransverse-Colon2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Sigmoid-Colon_Sigmoid-Colon_HCLAdultSigmoid-Colon1/Sigmoid-Colon_Sigmoid-Colon_HCLAdultSigmoid-Colon1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine1/Intestine_Intestine_HCLFetalIntestine1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine2/Intestine_Intestine_HCLFetalIntestine2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine3/Intestine_Intestine_HCLFetalIntestine3.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine4/Intestine_Intestine_HCLFetalIntestine4.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Intestine_Intestine_HCLFetalIntestine5/Intestine_Intestine_HCLFetalIntestine5.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Rectum_Rectum_HCLAdultRectum1/Rectum_Rectum_HCLAdultRectum1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Bladder_PartialTissue_HCLAdultBladder1/Bladder_PartialTissue_HCLAdultBladder1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Bladder_PartialTissue_HCLAdultBladder2/Bladder_PartialTissue_HCLAdultBladder2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Ureter_LowerUreter_HCLAdultUreter/Ureter_LowerUreter_HCLAdultUreter.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Duodenum_Duodenum_HCLAdultDuodenum1/Duodenum_Duodenum_HCLAdultDuodenum1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Ileum_Ileum_HCLAdultIleum2/Ileum_Ileum_HCLAdultIleum2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Jejunum_Jejunum_HCLAdultJejunum2/Jejunum_Jejunum_HCLAdultJejunum2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Kidney_AdjacentNormalTissue_HCLAdultKidney2/Kidney_AdjacentNormalTissue_HCLAdultKidney2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Kidney_PartialTissue_HCLAdultKidney3/Kidney_PartialTissue_HCLAdultKidney3.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Bronchus_BronchialEpithelialCell_Plasschaert2018/Bronchus_BronchialEpithelialCell_Plasschaert.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/heart_heart_cui2019/heart_heart_cui2019.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_Braga2019/lung_lung_Braga2019.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_Madissoon2019/lung_lung_Madissoon2019.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_carcinomas2019/lung_lung_carcinomas2019.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/heart_heart_asp2019cellpress/heart_heart_asp2019cellpress.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/heart_heart_HCLAdultHeart1/heart_heart_HCLAdultHeart1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/heart_heart_HCLAdultHeart2/heart_heart_HCLAdultHeart2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/heart_heart_HCLFetalHeart1/heart_heart_HCLFetalHeart1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/heart_heart_HCLFetalHeart2/heart_heart_HCLFetalHeart2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Heart_Heart_Wang2020/Heart_Heart_Wang2020.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Heart_Heart_Wang2020_2/Heart_Heart_Wang2020_2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_HCLAdultLung1/lung_lung_HCLAdultLung1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_HCLAdultLung2/lung_lung_HCLAdultLung2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_HCLAdultLung3/lung_lung_HCLAdultLung3.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_HCLFetalLung1/lung_lung_HCLFetalLung1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/lung_lung_HCLFetalLung2/lung_lung_HCLFetalLung2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Blood_Blood_HCLCord-Blood1/Blood_Blood_HCLCord-Blood1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Blood_Blood_HCLCord-Blood2/Blood_Blood_HCLCord-Blood2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood1/Blood_Blood_HCLAdult-Peripheral-Blood1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood2/Blood_Blood_HCLAdult-Peripheral-Blood2.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood3/Blood_Blood_HCLAdult-Peripheral-Blood3.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Blood_Blood_HCLAdult-Peripheral-Blood4/Blood_Blood_HCLAdult-Peripheral-Blood4.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow1/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/Spleen_Spleen_HCLAdult-Spleen1/Spleen_Spleen_HCLAdult-Spleen1.seuratobj.dbupload_v1.rds",
#   "/stor/public/hcad/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow2/BoneMarrow_BoneMarrow_HCLAdult-Bone-Marrow2.seuratobj.dbupload_v1.rds"
# )