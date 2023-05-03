RunAnchorIntegration=function(object_list,npcs=10,features=NULL){
  # select features that are repeatedly variable across datasets for integration
  if(is.null(features)){
	features <- SelectIntegrationFeatures(object.list = object_list)
  }
  anchors <- FindIntegrationAnchors(object.list = object_list, anchor.features = features)
  # this command creates an 'integrated' data assay
  combined <- IntegrateData(anchorset = anchors)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs =npcs, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:npcs)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:npcs)
  combined <- FindClusters(combined, resolution = 0.5)
  return(combined)
}

