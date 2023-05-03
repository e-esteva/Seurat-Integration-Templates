SeuratHarmonyIntegration=function(object_list,vars){
  obj.dims=lapply(object_list,function(x) dim(x)[2])
  
  object_list_counts=do.call('cbind',lapply(object_list,function(x) x@assays$RNA@counts))
  
  integrated_ <- CreateSeuratObject(counts = object_list_counts, project = "integrated_", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = integrated_@var.genes, npcs = 20, verbose = FALSE)
  
  integrated_@meta.data$dataset <- unlist(lapply(seq(length(vars)),function(x) rep(vars[x],obj.dims[x])))
  
  integrated_ <- integrated_ %>% 
    RunHarmony('dataset', plot_convergence = F)
  
  return(integrated_)
  
}
