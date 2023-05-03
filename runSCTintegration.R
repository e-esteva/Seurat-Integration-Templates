RunSCT_integration=function(query,ref,npcs=20,nfeatures=3e3){
  
  ref=SCTransform(ref, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE) 
  
  query=SCTransform(query, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE) 
  
  int.list <- list(ref = ref, query = query)
  features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = nfeatures)
  int.list <- PrepSCTIntegration(object.list = int.list, anchor.features = features)
  
  integration.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT",
                                           anchor.features = features)
  integration.combined.sct <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT")
  return(integration.combined.sct)
  
}