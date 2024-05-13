## Function to perform ORA (see shiny_run)

run_ora<-function(dataset.name, db.name, output.name="run"){
  set.seed(1234)
  # file.prefix and output dir
  file.prefix <- strsplit(dataset.name,"\\.")[[1]][1] #remove ext if there
  output.dir <- file.path("../",output.name, file.prefix)
  
  # Retrieve params
  par.fn <- paste0(file.prefix, "__ora_params.rds")
  params <- readRDS(file.path(output.dir, "ora",par.fn))
  org.db.name <- params$org.db.name
  print(org.db.name)
  minGSSize <- params$minGSSize
  maxGSSize <- params$maxGSSize

  
  # Retrieve geneList 
  gl.fn <- paste0(file.prefix, "__ora_input.rds")
  geneList <- readRDS(file.path(output.dir, "ora",gl.fn))
  
  geneKEGG <- geneList %>%
    dplyr::filter(ora.set == 1) %>%
    dplyr::pull(ENTREZID)
  universeKEGG <- pull(geneList, ENTREZID)
  universeKEGG <- as.character(universeKEGG) #not integers


  geneGO <- geneList %>%
    dplyr::filter(ora.set == 1) %>%
    dplyr::pull(SYMBOL)
  universeGO <- pull(geneList, SYMBOL)
  universeGO <- as.character(universeGO) #not integers

  # Use entire genome if not a susbset (i.e., by size or by missing p.adjvalue)
  if(!'p.adjvalue' %in% names(geneList) | !length(universeKEGG) > length(geneKEGG)){
    universeKEGG <- NULL
  }
  if(!'p.adjvalue' %in% names(geneList) | !length(universeGO) > length(geneGO)){
    universeGO <- NULL
  }


  if(db.name=='db_kegg'){
    # Perform ORA
  enrichResult <- clusterProfiler::enrichKEGG(
    geneKEGG,
    organism = "mmu",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    keyType = 'ncbi-geneid',
    # pAdjustMethod="holm", #default is "BH"
    pvalueCutoff = 0.05 #to limit results
    )
  }
  else if(db.name=='db_go_bp') {
enrichResult <- enrichGO(
                           geneGO,
                           OrgDb = org.db.name,
                           keyType = 'SYMBOL',
                           readable = T,
                           ont = "BP",
                           pvalueCutoff = 0.05)

  }
    else if(db.name=='db_go_mf') {
enrichResult <- enrichGO(
                           geneGO,
                           OrgDb = org.db.name,
                           keyType = 'SYMBOL',
                           readable = T,
                           ont = "MF",
                           pvalueCutoff = 0.05)

  }
    else if(db.name=='db_go_cc') {
enrichResult <- enrichGO(
                           geneGO,
                           OrgDb = org.db.name,
                           keyType = 'SYMBOL',
                           readable = T,
                           ont = "CC",
                           pvalueCutoff = 0.05)

  }
  else{
  # Perform ORA
    # Object from string
  database <- eval(parse(text=db.name))
  enrichResult <- clusterProfiler::enricher(
    geneKEGG,
    universe = universe,
    TERM2GENE = database[,c("term","gene")],
    TERM2NAME = database[,c("term","name")],
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    # pAdjustMethod="holm", #default is "BH"
    pvalueCutoff = 0.05 #to limit results
    )
    }
  if(!is.null(enrichResult))
    enrichResult <- setReadable(enrichResult, eval(parse(text=org.db.name)), keyType = "ENTREZID")
  
  #Prune for size
  enrichResult@geneSets <- list()
  
  # Save objects
  gl.fn <- paste(file.prefix, db.name,"ora","geneList.rds", sep = "_")
  saveRDS(geneKEGG, file.path(output.dir,"ora",gl.fn))
  er.fn <- paste(file.prefix, db.name,"ora","result.rds", sep = "_")
  saveRDS(enrichResult, file.path(output.dir,"ora",er.fn))
  
  ## Save df as TSV and XLSX
  enrichResult.df <- as.data.frame(enrichResult)
  tsv.fn <- paste(file.prefix, db.name,"ora.tsv", sep = "_")
  xlsx.fn <- paste(file.prefix, db.name,"ora.xlsx", sep = "_")
  write.table(enrichResult.df,file.path(output.dir,"ora",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(enrichResult.df,file.path(output.dir,"ora",xlsx.fn))
  
  ## Plot
  if (nrow(enrichResult.df) > 0)
    plot_results(enrichResult, geneKEGG, file.prefix, output.dir, db.name, "ora")

}
