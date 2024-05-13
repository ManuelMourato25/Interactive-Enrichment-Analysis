## Function to perform GSEA (see shiny_run)

run_gsea<-function(dataset.name, db.name, output.name="run"){
  set.seed(1234)
  # file.prefix and output dir
  file.prefix <- strsplit(dataset.name,"\\.")[[1]][1] #remove ext if there
  output.dir <- file.path("../",output.name, file.prefix)
  
  # Retrieve params
  par.fn <- paste0(file.prefix, "__gsea_params.rds")
  params <- readRDS(file.path(output.dir, "gsea",par.fn))
  org.db.name <- params$org.db.name
  minGSSize <- params$minGSSize
  maxGSSize <- params$maxGSSize
  

  
  # geneList from file.prefix
  gl.fn <- paste0(file.prefix, "__gsea_input.rds")
  geneList <- readRDS(file.path(output.dir, "gsea",gl.fn))

  # Sorted named list for clusterProfiler, a.k.a. geneList
  ## pre-ranking
  # geneList <- geneList %>%
  #   mutate(rank = rank(rank,  ties.method = "random")) %>%
  #   arrange(desc(rank))
  ranked.genes.entrez.nl<-geneList$rank
  names(ranked.genes.entrez.nl)<-geneList$ENTREZID
  geneListKEGG <- sort(ranked.genes.entrez.nl, decreasing = T)

  ranked.genes.symbol.nl<-geneList$rank
  names(ranked.genes.symbol.nl)<-geneList$SYMBOL
  geneListGO <- sort(ranked.genes.symbol.nl, decreasing = T)
  
  # Perform GSEA
  if (db.name=='db_kegg'){

       gseaResult <- gseKEGG(geneList= geneListKEGG,
               organism     = 'mmu',
               nPerm        = 10000,
               minGSSize    = minGSSize,
               maxGSSize    = maxGSSize,
               pvalueCutoff = 0.05,
               keyType       = "ncbi-geneid")
       if(!is.null(gseaResult))
          gseaResult <- setReadable(gseaResult, eval(parse(text=org.db.name)), keyType = "ENTREZID")
     }
       else if (db.name=='db_go_bp'){

       gseaResult <- gseGO(geneList=geneListGO,
                     ont ="BP",
                     keyType = "SYMBOL",
                     nPerm = 10000,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = org.db.name,
                     pAdjustMethod = "BH")
     }
           else if (db.name=='db_go_mf'){

       gseaResult <- gseGO(geneList=geneListGO,
                     ont ="MF",
                     keyType = "SYMBOL",
                     nPerm = 10000,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = org.db.name,
                     pAdjustMethod = "BH")
     }
           else if (db.name=='db_go_cc'){

       gseaResult <- gseGO(geneList=geneList,
                     ont ="CC",
                     keyType = "SYMBOL",
                     nPerm = 10000,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = org.db.name,
                     pAdjustMethod = "BH")
     }
     else{
            # Object from string
              database <- eval(parse(text=db.name))
              gseaResult <- clusterProfiler::GSEA(
                geneList,
                TERM2GENE = database[,c("term","gene")],
                TERM2NAME = database[,c("term","name")],
                minGSSize = minGSSize,
                maxGSSize = maxGSSize,
                # pAdjustMethod="holm", #default is "BH"
                pvalueCutoff = 0.05, #to limit results
                verbose=FALSE)
     }
  
  # Save objects
  gl.fn <- paste(file.prefix, db.name,"gsea","geneList.rds", sep = "_")
  saveRDS(geneList, file.path(output.dir,"gsea",gl.fn))
  er.fn <- paste(file.prefix, db.name,"gsea","result.rds", sep = "_")
  saveRDS(gseaResult, file.path(output.dir,"gsea",er.fn))
  
  ## Save df as TSV and XLSX
  gseaResult.df <- as.data.frame(gseaResult)
  tsv.fn <- paste(file.prefix, db.name,"gsea.tsv", sep = "_")
  xlsx.fn <- paste(file.prefix, db.name,"gsea.xlsx", sep = "_")
  write.table(gseaResult.df,file.path(output.dir,"gsea",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(gseaResult.df,file.path(output.dir,"gsea",xlsx.fn))
  
  ## Plot
  if (nrow(gseaResult.df  ) > 0)
    plot_results(gseaResult, geneList, file.prefix, output.dir, db.name, "gsea")
}
